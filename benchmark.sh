#! /bin/bash

julia --project=. --startup-file=no -e 'import Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=. --startup-file=no -e 'using BenchmarkTools; include("src/exact.jl"); @btime exact()'
julia --project=. --startup-file=no -e 'using BenchmarkTools; include("src/fdm.jl"); @btime fdm()'
