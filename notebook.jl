### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ cf07621a-9a65-11ec-1fc0-51fbadbb086f
using Pkg; Pkg.activate()

# ╔═╡ 345a90f7-f8ac-43b2-b5bd-0282d9447556
begin
	using DifferentialEquations
	using ModelingToolkit
	using DynamicalSystems
	using ForwardDiff
	using LinearAlgebra
	using Graphs
	using Plots
	using PlutoUI
	using LaTeXStrings
end

# ╔═╡ eaf9dda2-d898-4771-9f86-ca4be174712c
md"""
# Dynamical Systems in Julia
In this tutorial, we'll looking at specifying, simulating and analyzing dynamical systems in Julia. This is one area where Julia especially shines, thanks in no small part to the `DifferentialEquations` package, which pretty much provides the largest and fastest suite of ODE solvers on the market today.
"""

# ╔═╡ eb6aed5a-0f95-4a76-8a99-b8719984bbe5
md"## Defining systems symbolically"

# ╔═╡ 0e28ede5-0af6-49a7-b95d-7094e320eee4
md"""
We'll start by using ModelingToolkit to define the famous Lorenz system. The process starts by a symbolic description of the problem.
"""

# ╔═╡ 6ccfc87f-40e7-4ee4-baae-ee41b81efe67
begin
	@parameters t σ ρ β
	@variables x(t) y(t) z(t)
	D = Differential(t)

	eqs = [
		D(x) ~ σ*(y-x),
       	D(y) ~ x*(ρ-z)-y,
       	D(z) ~ x*y - β*z
	]

	@named lorenz = ODESystem(eqs)
end

# ╔═╡ 3e4ec586-482f-4c1b-a8e7-4ae88e48fa38
md"""
So far, the `lorenz` is purely symbolic, which has the added benefit of being able to render the equations in LaTeX as above! The next step before simulating it is defining values for our initial condition and parameters.
"""

# ╔═╡ ce0b65b5-e57a-4def-acb6-892c49f182d9
begin
	u0 = [
		x => 0.0,
      	y => 0.5,
      	z => 0.5
	]

	p  = [
		ρ => 28.0,
      	σ => 10.0,
      	β => 8/3
	]
end

# ╔═╡ 1a7e8273-8069-42a8-ab7b-62f5fe1ade67
md"""
Simulating our system is a two step process. After defining our timespan of integration, we first create an `ODEProblem`, we describes an ODE + initial conditions and parameters, and then pass it to the `solve` function to compute the solution of our problem. The second argument to solve is the solver we want to use (`Tsit5` is a 4-5 order Runge-Kutta scheme which is generally a good default choice for non-stiff problems)
"""

# ╔═╡ 1d77a2be-773c-4edd-9d68-1ed84fb5dadf
begin
	tspan = (0.0,100.0)
	prob = ODEProblem(lorenz,u0,tspan,p,jac=true)
	sol = solve(prob,Tsit5())
end

# ╔═╡ 7c712f73-3931-4cf4-82d8-9d9e2263df14
plot(sol, vars=(x,y,z))

# ╔═╡ 432fa70e-9225-4e76-aec3-c27d31887fb9
md"""
NB. There are also packages for modeling specific classes of dynamical systems, such as `Catalyst` for chemical reaction networks, or `PowerDynamics` for power grids.
"""

# ╔═╡ 60063cc3-e71f-461b-be0a-491b1b566718
md"## Defining Systems by hand"

# ╔═╡ ea2eeac2-5514-4fe8-8921-6235cdb0d1ad
md"""
While defining systems symbolically is cool, this can be tedious for large systems, and sometimes you just want full control over how your ODE function is defined. To that end, we will now define another famous system, the Brusselator by hand and then build a network reaction diffusion system with it.

The Brusselator equations are as follows,
```math
\begin{align*}
\dot{x} &= a + x^2y -bx -x \\
\dot{y} &= bx - x^2y
\end{align*}
```
"""

# ╔═╡ 2bf13604-71c6-474f-90de-5b270c5c0708
md"""
We can construct ODEProblems by providing a function that computes the derivative, just like we do in Matlab. The functions we can use come in two flavours. The first type just takes as input the state `u`, parameters stored in `p` and time `t` (for non autonomous systems), and returns the new state.
"""

# ╔═╡ 743b31cd-6b2a-472c-8a70-a3955a2c26b6
function brusselator(u, p, t)
	a, b = p
	[
		a + u[1]^2 * u[2] - b*u[1] - u[1],
		b*u[1] - u[1]^2 * u[2]
	]
end

# ╔═╡ 03bd6ff6-107f-4b04-9cb8-fe5c28570e58
md"""
However, it is always recommended to use the second form, which also takes as input a preallocated array for the output and modifies it. Doing this is more efficient in terms of memory allocation, which is actually an important source of slowdowns.
"""

# ╔═╡ 08fc4616-d2af-471b-adc0-195bb566334f
function brusselator!(du, u, p, t)
	a, b = p
	du[1] = a + u[1]^2 * u[2] - b*u[1] - u[1]
	du[2] = b*u[1] - u[1]^2 * u[2]
end

# ╔═╡ d0c1e91e-d1df-4a12-a0e2-8513eeee7ba6
begin
	u0_b = [1.0, 1.0]
	p_b = [1.0, 3.0]
end

# ╔═╡ 8623c9b7-5724-4678-9b53-4312389a276d
begin
	prob_b = ODEProblem(brusselator!, u0_b, (0.0, 50.0), p_b)
	sol_b = solve(prob_b, Tsit5())
end

# ╔═╡ da0220ea-99aa-4a0d-8826-f0e461fd503e
plot(sol_b)

# ╔═╡ 60d1ea9b-b935-4e88-9e92-07497f06aeaa
md"""
Now that we have our brusselator function, let's use it on a network to make a reaction diffusion system. Such a system is of the form

```math
\dot{u_i} = f(u_i) + D \sum_{j=1}^N A_{ij}(u_j - u_i),
```
where ``u_i\in\mathbb{R}^m`` is the state of node ``i``, subject to local dynamics ``f:\mathbb{R}^m\rightarrow \mathbb{R}^m``, ``A\in\mathbb{R}^{N\times N}`` is the adjacency matrix of the graph (i.e. ``A_{ij}=1`` iff there is a link between nodes ``i`` and ``j``), and ``D`` is a diagonal matrix of diffusion coefficients.

We can rewrite this system in vector form (glossing over some details) as

```math
\dot{\mathbf{u}} = F(\mathbf{u}) + \mathbf{D} L \mathbf{u},
```
where ``L`` is the *Laplacian matrix* of the graph.

Reaction-Diffusion systems can take more general forms, but we'll stick to this one for today.
"""

# ╔═╡ efbb46be-c442-4ab3-85a0-b10cbc3a2a9e
md"""
Now, we need to create a function for our reaction-diffusion system.
"""

# ╔═╡ 8b2ec98a-b37f-4bd0-9fc8-9f3dcf9550c9
function brusselator_rd(u, p, t)
	a, b, L, D = p
	du = similar(u) # creates an array of same shape and element type as u
	for i in 1:size(du, 1)
		du[i,:] .= brusselator(u[i,:], [a, b], t)
	end
	du[:,1] .-= D[1] .* (L*u[:,1])
	du[:,2] .-= D[2] .* (L*u[:,2])
	return du
end

# ╔═╡ 25c22a1f-6ae4-4b72-8e4c-de844e2699d1
function brusselator_rd!(du, u, p, t)
	# NB. the commented out parts should be equivalent to what is used, but for some reason it is not
	a, b, L, D = p
	du .= 0.0
	for i in 1:size(du, 1)
		#brusselator!(du[i,:], u[i,:], [a, b], t)
		du[i,:] .+= brusselator(u[i,:], [a, b], t)
	end
	#mul!(du[:,1], L, u[:,1], -D[1], 0.0) # du[:,1] .+= D[1]*L*u[:,1]
	#mul!(du[:,2], L, u[:,2], -D[2], 0.0)
	du[:,1] .-= D[1] .* (L*u[:,1])
	du[:,2] .-= D[2] .* (L*u[:,2])
end

# ╔═╡ 582b74c7-aaf6-4a8e-901f-489199ccc0d6
md"""
Now, we need a network to run our system on. To do this, we can use the `Graphs` package.
"""

# ╔═╡ 3e575a45-50d6-4297-84e5-6699e36d00dd
g = erdos_renyi(10, 0.4)

# ╔═╡ ed1de34a-a42e-474d-ba36-411163a343f8
md"""
We can use the function `laplacian_matrix` to get the Laplacian matrix of the graph. The nice thing about this function is that it returns a sparse matrix, which uses optimized methods for matrix multiplications.
"""

# ╔═╡ 947c169c-a900-4728-ad29-5200b52f9599
laplacian_matrix(g)

# ╔═╡ 04bbc594-17c9-4b41-a599-d41e2e541c9f
begin
	u_rd_0 = 5*rand(nv(g), 2)
end

# ╔═╡ 29a0f8c8-5529-496c-99f1-b436d15f2517
md"""
Here's a nice feature of Pluto. It provides allows interactively modifying variables with html elements (sliders, buttons, ...) via the `PlutoUI` package.

Let's use it to make sliders for the parameters of our system.
"""

# ╔═╡ f346a0b3-4d1d-49ad-94ec-cc42ee75a05b
md"""
a: $(@bind a Slider(0.0:0.01:5.0, default=1.0, show_value=true))

b: $(@bind b Slider(0.0:0.01:5.0, default=3.0, show_value=true))

D₁: $(@bind D_1 Slider(0.0:0.01:10.0, default=1.0, show_value=true))

D₂ $(@bind D_2 Slider(0.0:0.01:10.0, default=1.0, show_value=true))
"""

# ╔═╡ 8a1c55e5-1008-44ee-b853-c53c2dbf3bb2
p_rd = [a, b, laplacian_matrix(g), [D_1, D_2]];

# ╔═╡ 2ee789c5-ae38-4e85-95e1-920a580df058
begin
	prob_rd = ODEProblem(brusselator_rd!, u_rd_0, (0.0, 50.0), p_rd)
	sol_rd = solve(prob_rd, Tsit5())
	nothing
end

# ╔═╡ 9d1f3408-e0f5-4093-bb6b-862485e99476
begin
	local p1 = plot(sol_rd, vars=(0,1:nv(g)), label="", ylabel=L"x_i")
	local p2 = plot(sol_rd, vars=(0,(nv(g)+1:2nv(g))), label="", ylabel=L"y_i")
	plot(p1, p2)
end

# ╔═╡ a4d4a26f-ff25-4bf0-b8ab-2db56131ee3e
md"""
### Exercise
Implement your favourite dynamical system and play around with it's parameters.
"""

# ╔═╡ 91f1ce37-0d68-4ed2-9359-dbcfd06f6cdd
md"""
## Intermezzo: using `DynamicalSystems` to study chaotic Systems

Before moving to the last part of this tutorial, let us make a detour to look at the `DynamicalSystems` package, which provides a library of famous dynamical systems, as well as generic tools for ananlyzing them, particularly chaotic systems. 
"""

# ╔═╡ 4328e814-82a8-4e7d-8b8a-44a9d5ba17c5
ds = Systems.logistic()

# ╔═╡ 7b11ac9b-a470-4059-bbad-ba6bc443c68b
md"""
Here we load the well known logistic map from the library of dynamical systems, and generate its orbit diagram for varying values of the `r` parameter. The magic function is `orbitdiagram`, which should work for pretty much any discrete map.
"""

# ╔═╡ a5629c28-53bf-4b19-a49f-ab26df1d1b46
begin
	i = 1
	pvalues = 1:0.005:4
    ics = [rand() for m in 1:10]
    n = 2000
    Ttr = 2000
    p_index = 1
	
    output = orbitdiagram(ds, i, p_index, pvalues; n = n, Ttr = Ttr)
end

# ╔═╡ 74adb98a-1ad5-4181-a3be-ec5fe2565a35
begin
	local L = length(pvalues)
	local x = Vector{Float64}(undef, n*L)
	local y = copy(x)
	for j in 1:L
    	x[(1 + (j-1)*n):j*n] .= pvalues[j]
  		y[(1 + (j-1)*n):j*n] .= output[j]
	end
	scatter(x,y,
		xlabel="r",
		ylabel="x",
		label="",
		markersize=0.8,
		color=:black,
		markeralpha=0.05,
		markerstrokealpha=0.0
	)
end

# ╔═╡ c38e202f-a2d2-4449-8308-df22774b7af1
md"""
NB. While `DynamicalSystems` provides a bunch of cool tools, a lot of them seem to inexplicably break when used outside of the examples. For example, try changing the above code to use values of `r` above 4.
"""

# ╔═╡ d5ab9847-ace4-4a5b-8c81-d2eb8b56c896
md"""
## Analyzing the stability of our system

Reaction-Diffusion systems like the one we just saw are interesting because of their asymptotic behaviour.  Because the laplacian always has a zero eigenvalue associated to an eigenvector whose entries are all the same, a reaction-diffusion system always admits a homogeneous solution where the nodes are all synchronised.
By tuning the diffusion coefficients, it's possible to (de)synchronize the states of the nodes.

The standard tool for analysing this phenomenon is the Master Stability Function [Pecora-Caroll. 1998].

The main idea consists in expressing the system in the basis of eigenvectors of the laplacian matrix. By considering a small perturbation around the synchronized solution, and expressing its linearized equations in the eigenbasis of ``L``, one can obtain a set of decoupled equations involving the laplacian eigenvalues and the jacobian of the local dynamics.

Skipping over the details, we can check that the synchronized solution is linearly stable by checking that the following ODE is stable for every eigenvalue ``\lambda`` of the laplacian matrix

```math
\dot{\xi} = [\mathcal{J}f(u_*) - \lambda D] \xi
```
where ``\mathcal{J}f`` denotes the jacobian of ``f``.

This can be done, for example, by computing the Lyapunov exponent of this equation, thereby yielding a function ``MSF(\lambda)`` called the *Master Stability Function* of the system.
"""

# ╔═╡ 575e3106-704e-4bea-8155-7ee6e6c9da5d
md"First, let's compute the eigenvalues of the laplacian matrix. We don't actually need to compute the eigenvectors."

# ╔═╡ bff34b78-f7ad-452e-9b95-a498f3e377f5
λ = laplacian_spectrum(g)

# ╔═╡ bba5eaa1-16f2-401d-8cf3-f05e9761583d
md"""
Next, let's start by writing a function that takes the local dynamics as input, and generates a function for the MSF ODE.

Notice the use, of `ForwardDiff.jacobian` to automagically compute the Jacobian matrix using Automatic Differentiation (this is different from Finite Differences or Symbolic methods)
"""

# ╔═╡ ddbb84d1-945a-4184-850b-31ee3540e960
function msf_ode_fun(f)
	function msf_ode!(du, u, p, t)
		p_f, D, λ = p
		du[:,1] .= f(u[:,1], p_f, t)
		J = ForwardDiff.jacobian(x -> f(x,p_f,t), u[:,1])
		du[:,2] .= J*u[:,2] - λ*D*u[:,2]
	end
	return msf_ode!
end

# ╔═╡ f0215a78-146d-4009-ae73-8a2f35072db9
md"""
Now, we can move on the computing the MSF proper. To do this, we solve the MSF ODE, then estimating the Lyapunov exponent from the time series.

NB. Recall that the Lyapunov exponent for a dynamical sytem is defined by
```math
\lambda = \lim_{t\rightarrow \infty} \lim_{\xi_0 \rightarrow 0} \frac{1}{t} \ln \frac{\|\xi(t)\|}{\|\xi_0\|},
```
where ``\xi(0) = \xi_0`` is a small perturbation applied to the orbit of interest (here, the perdiodic solution).
"""

# ╔═╡ 84d86795-e2c4-4c2c-888b-ae999bcda2dd
function master_stability_function(f, u_0, p_f, D, t; δ=0.1, k=50, Δt=0.1)
	msf_ode! = msf_ode_fun(f)
	n = length(u_0)
	u_msf_0 = hcat(u_0, δ*randn(n))
	function msf(λ; u_msf_0=u_msf_0, p_f=p_f, D=D, t=t, Δt=Δt)
		p_msf = [p_f, D, λ]
		prob = ODEProblem(msf_ode!, u_msf_0, (0.0,t), p_msf)
		sol = solve(prob, Tsit5(), saveat=0.0:Δt:t, abstol=1e-10, reltol=1e-10)
		ξ_0 = norm(sol.u[1][:,2])
		l = length(sol.u)
		ξ_t = norm([sol.u[k][:,2] for k in l-k:l])
		return log(ξ_t/ξ_0)/sol.t[end]
	end
	return msf
end

# ╔═╡ 8c592351-4da7-4bed-a99c-595156c1ed8b
begin
	u0_msf = sol_b.u[end]
	p_f = [1.0, 3.0]
	D_msf = Diagonal([0.2,5.0])
end

# ╔═╡ e3411cdb-acf7-40da-9f5e-6852c1370797
msf = master_stability_function(brusselator, u0_msf, p_f, D_msf, 200.0)

# ╔═╡ 5fbfe4aa-f2e5-4845-92e7-9de00063059f
msf(0.0)

# ╔═╡ ffecadd1-95b5-412f-86e4-e4dee733cab1
md"""
Note that this implementation is not super fast, as it takes about a minute to run the next cell. There's lots of room for improvements though, and in the mean time, we still have a generic funcion for computing the MSF.
"""

# ╔═╡ bad93338-72e3-46fc-888f-e507c62cefc6
begin
	local λs = 0.0:0.01:7.0
	plot(λs, msf.(λs), 
		label="",
		xlabel="-Re λ"
	)
	#plot!(λs, msf_ds.(λs), label="ds")
	plot!([0.0,7.0],[0.0,0.0], linestyle=:dash, label="", color=:black)
	scatter!(λ, msf.(λ), label="")
end

# ╔═╡ 63c14c9e-808e-464b-b5cc-25921aa34f56
md"""
Below, we examine the trajectory of the MSF equations for λ = 0. Note that this is quite ill-behaved (notice the log scale), so it's no wonder this is tricky to get right.
"""

# ╔═╡ 8f8a9dec-4de2-4fc7-9d87-a3638ec47ee3
begin
	f! = msf_ode_fun(brusselator)
	u0_ = hcat(u0_msf, 0.1*randn(2))
	p_msf = [p_f, D_msf, 0.0]
	prob_msf = ODEProblem(f!, u0_, (0.0,100.0), p_msf)
	sol_msf = solve(prob_msf, Rodas5(), abstol=1e-10, reltol=1e-10)
end

# ╔═╡ e85b9b6d-736d-446a-acb1-b6ed0ccdd907
begin
	local foo = (t,u,v) -> (t,norm([u;v]))
	plot(sol_msf, vars=(foo,0,3,4), 
		yscale=:log10, 
		label="", 
		ylabel=L"\ln \|\xi\|"
	)
end

# ╔═╡ Cell order:
# ╟─eaf9dda2-d898-4771-9f86-ca4be174712c
# ╠═cf07621a-9a65-11ec-1fc0-51fbadbb086f
# ╠═345a90f7-f8ac-43b2-b5bd-0282d9447556
# ╟─eb6aed5a-0f95-4a76-8a99-b8719984bbe5
# ╟─0e28ede5-0af6-49a7-b95d-7094e320eee4
# ╠═6ccfc87f-40e7-4ee4-baae-ee41b81efe67
# ╟─3e4ec586-482f-4c1b-a8e7-4ae88e48fa38
# ╠═ce0b65b5-e57a-4def-acb6-892c49f182d9
# ╟─1a7e8273-8069-42a8-ab7b-62f5fe1ade67
# ╠═1d77a2be-773c-4edd-9d68-1ed84fb5dadf
# ╠═7c712f73-3931-4cf4-82d8-9d9e2263df14
# ╟─432fa70e-9225-4e76-aec3-c27d31887fb9
# ╟─60063cc3-e71f-461b-be0a-491b1b566718
# ╟─ea2eeac2-5514-4fe8-8921-6235cdb0d1ad
# ╟─2bf13604-71c6-474f-90de-5b270c5c0708
# ╠═743b31cd-6b2a-472c-8a70-a3955a2c26b6
# ╟─03bd6ff6-107f-4b04-9cb8-fe5c28570e58
# ╠═08fc4616-d2af-471b-adc0-195bb566334f
# ╠═d0c1e91e-d1df-4a12-a0e2-8513eeee7ba6
# ╠═8623c9b7-5724-4678-9b53-4312389a276d
# ╠═da0220ea-99aa-4a0d-8826-f0e461fd503e
# ╟─60d1ea9b-b935-4e88-9e92-07497f06aeaa
# ╟─efbb46be-c442-4ab3-85a0-b10cbc3a2a9e
# ╠═8b2ec98a-b37f-4bd0-9fc8-9f3dcf9550c9
# ╠═25c22a1f-6ae4-4b72-8e4c-de844e2699d1
# ╟─582b74c7-aaf6-4a8e-901f-489199ccc0d6
# ╠═3e575a45-50d6-4297-84e5-6699e36d00dd
# ╟─ed1de34a-a42e-474d-ba36-411163a343f8
# ╠═947c169c-a900-4728-ad29-5200b52f9599
# ╠═04bbc594-17c9-4b41-a599-d41e2e541c9f
# ╟─29a0f8c8-5529-496c-99f1-b436d15f2517
# ╟─f346a0b3-4d1d-49ad-94ec-cc42ee75a05b
# ╠═8a1c55e5-1008-44ee-b853-c53c2dbf3bb2
# ╠═2ee789c5-ae38-4e85-95e1-920a580df058
# ╠═9d1f3408-e0f5-4093-bb6b-862485e99476
# ╟─a4d4a26f-ff25-4bf0-b8ab-2db56131ee3e
# ╟─91f1ce37-0d68-4ed2-9359-dbcfd06f6cdd
# ╠═4328e814-82a8-4e7d-8b8a-44a9d5ba17c5
# ╟─7b11ac9b-a470-4059-bbad-ba6bc443c68b
# ╠═a5629c28-53bf-4b19-a49f-ab26df1d1b46
# ╠═74adb98a-1ad5-4181-a3be-ec5fe2565a35
# ╟─c38e202f-a2d2-4449-8308-df22774b7af1
# ╟─d5ab9847-ace4-4a5b-8c81-d2eb8b56c896
# ╟─575e3106-704e-4bea-8155-7ee6e6c9da5d
# ╠═bff34b78-f7ad-452e-9b95-a498f3e377f5
# ╟─bba5eaa1-16f2-401d-8cf3-f05e9761583d
# ╠═ddbb84d1-945a-4184-850b-31ee3540e960
# ╟─f0215a78-146d-4009-ae73-8a2f35072db9
# ╠═84d86795-e2c4-4c2c-888b-ae999bcda2dd
# ╠═8c592351-4da7-4bed-a99c-595156c1ed8b
# ╠═e3411cdb-acf7-40da-9f5e-6852c1370797
# ╠═5fbfe4aa-f2e5-4845-92e7-9de00063059f
# ╟─ffecadd1-95b5-412f-86e4-e4dee733cab1
# ╠═bad93338-72e3-46fc-888f-e507c62cefc6
# ╟─63c14c9e-808e-464b-b5cc-25921aa34f56
# ╠═8f8a9dec-4de2-4fc7-9d87-a3638ec47ee3
# ╠═e85b9b6d-736d-446a-acb1-b6ed0ccdd907
