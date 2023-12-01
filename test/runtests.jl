using ADAPT
using Test

@testset "ADAPT.jl" begin
    # TEMP: this are not proper tests. They are scripts. For now, pass if they run.
    include("hubbard_qeb.jl")
    include("check_unitary.jl")
end
