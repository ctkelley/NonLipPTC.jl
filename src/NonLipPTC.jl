module NonLipPTC

using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using SIAMFANLEquations.Examples
using PythonPlot
using LaTeXStrings
using Printf
import SIAMFANLEquations: Anderson_Init, EvalF!, Orthogonalize!,
    ItStatsA, updateHist!, falpha, updateStats!, AndersonOK,
    CloseIteration

include("Tools/fprintTeX.jl")
include("Tools/plothist.jl")
include("OneD.jl")
include("Alg_Test.jl")
include("Algorithms/Alg_Structs.jl")
include("Algorithms/alg1.jl")
include("Algorithms/alg3.jl")
include("Algorithms/alg1examples.jl")
include("test/fdgrad.jl")
include("Talk_Stuff/complexity.jl")

#include("Sanity.jl")
include("Figures/tau_test_copper.jl")
include("Figures/aa_test_copper.jl")
include("Problem_Data/boundary.jl")
include("Problem_Data/build_problemex1.jl")
include("Problem_Data/Exact.jl")
include("Problem_Data/NL_Equationex1.jl")

export alg1_test, alg3_test, alg0_test, alg0
export alg1e1, uefun_ex12d, tau_test_copper, aa_test_copper
export testP, aa_test, alg1e1, build_problem, FGen
export GF!, complexity
export fdgrad

end
