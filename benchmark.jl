using Printf
using BenchmarkTools

# Include source files

include("src/exact.jl")
include("src/fdm.jl")
include("src/qtt.jl")

# For precompilation

println("\n## Precompilation\n")
println("| Method | Energy |")
println("| ------ | ------ |")
println("| Exact   | ", exact(), " | ")
println("| FDM     | ", fdm(n = 2^10, x_min = -50.0, x_max = 50.0), " |")
println("| QTT+TCI | ", qtt(c =   10, x_min = -50.0, x_max = 50.0), " |")

# FDM v.s. QTT+TCI

println("\n## Energy\n")
println("|  c |        n |               FDM |           QTT+TCI |")
println("| -- | -------- | ----------------- | ----------------- |")
for c in 6:20
    n = 2^c
    @printf("| %2d | %8d | %.15f | %.15f |\n", c, n, fdm(n = n, x_min = -50.0, x_max = 50.0), qtt(c = c, x_min = -50.0, x_max = 50.0))
end

println("\n## Time\n")
println("|  c |        n |             FDM |         QTT+TCI |")
println("| -- | -------- | --------------- | --------------- |")
for c in 6:20
    n = 2^c
    fdm_bench = @benchmark fdm(n = $n, x_min = -50.0, x_max = 50.0)
    qtt_bench = @benchmark qtt(c = $c, x_min = -50.0, x_max = 50.0)
    # fdm_time = @sprintf("%3.1e ± %3.1e", mean(fdm_bench).time, std(fdm_bench).time)
    # qtt_time = @sprintf("%3.1e ± %3.1e", mean(qtt_bench).time, std(qtt_bench).time)
    @printf("| %2d | %8d | %12s | %12s |\n", c, n, BenchmarkTools.prettytime(mean(fdm_bench).time), BenchmarkTools.prettytime(mean(qtt_bench).time))
end