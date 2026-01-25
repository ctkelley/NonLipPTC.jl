module NonLipPTC

ENV["JULIA_CONDAPKG_VERBOSITY"] = "0"
ENV["CONDAPKG_QUIET"] = "1"
using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using SIAMFANLEquations.Examples
using PythonPlot;
using LaTeXStrings
using Printf

include("Tools/fprintTeX.jl")
include("Tools/plothist.jl")
include("OneD.jl")
include("Alg_Test.jl")
include("Algorithms/Alg_Structs.jl")
include("Algorithms/alg0.jl")
include("Algorithms/alg1.jl")
include("Algorithms/alg3.jl")

export OneD, POneD, setdata, proj, fobj1D
export alg1_test, alg3_test, alg0_test, alg0

end
