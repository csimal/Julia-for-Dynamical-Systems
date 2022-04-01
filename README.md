# Julia for Dynamical Systems
Resources for a tutorial about using Julia to Study Dynamical Systems given on 01/04/2022.

You may also want to check the page from my [previous talk](https://github.com/csimal/Julia-Unamur).

## Running the examples
The examples in the tutorial were run in a [Pluto notebook](https://github.com/fonsp/Pluto.jl). To run it, you'll need to install it, as well as a few packages. First install the latest stable release (v1.7.2) from [here](https://julialang.org/downloads/).

You can install the necessary packages easily by running the script `packages.jl`. Simply run
```
julia packages.jl
```
in a command line. (This will work on any OS, though on Windows you need to make sure to add Julia to the PATH during the installation).

Once you've done that (this can take some time), you can start a Pluto session by opening Julia in a command line, and running
```julia-repl
julia> using Pluto

julia> Pluto.run()
```
This will open a tab on your web browser where you can enter the path to where you downloaded `notebook.jl`. Once you've done that, the notebook will run itself.

## April Fool's joke
The notebook shown at the beginning is actually a legit tutorial on how to solve the Poisson equation using the package [Gridap.jl](https://gridap.github.io/Gridap.jl/stable/) which seems to provide a comprehensive toolkit for solving PDEs in Julia.

## DifferentialEquations
This package is single-handedly responsible for a good chunk of the growing adoption of Julia, as it basically provides the [fastest](https://github.com/SciML/SciMLBenchmarks.jl), and [most extensive](https://diffeq.sciml.ai/stable/) suite of ODE solvers available, making it *the* state of the art.

It stands as the core of the wider [SciML](https://sciml.ai/) package ecosystem which is all about using machine learning techniques for scientific research.

In the tutorial, we used the package [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) to define dynamical systems symbolically. There are also specific packages for modeling [chemical reactions](https://github.com/SciML/Catalyst.jl) and [power grids](https://juliaenergy.github.io/PowerDynamics.jl/stable/), and *many more*.

## DynamicalSystems
Simulating dynamical systems is not the whole story, we also want to analyse them. [DynamicalSystems](https://juliadynamics.github.io/DynamicalSystems.jl/latest/) is a toolbox of various methods for doing just that, especially for studying chaotic systems. It also doubles as a library of famous models, for more educational purposes. I've noticed that it tends to break when applied outside examples, so your mileage may vary. 
## Other Packages used in this tutorial

* [Graphs](https://github.com/JuliaGraphs/Graphs.jl)
* [Plots](https://github.com/JuliaPlots/Plots.jl)
* [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl)
* [LaTeXStrings](https://github.com/stevengj/LaTeXStrings.jl)



