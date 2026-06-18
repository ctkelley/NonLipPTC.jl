module NonLipPTC

using LinearAlgebra, SparseArrays
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
include("Tools/manefesto.jl")
include("OneD.jl")
#include("Alg_Test.jl")
include("Algorithms/Alg_Structs.jl")
include("Algorithms/alg1.jl")
include("Algorithms/alg2g.jl")
include("Algorithms/alg1examples.jl")
include("Algorithms/alg2gexamples.jl")
include("test/fdgrad.jl")
include("Talk_Stuff/complexity.jl")
include("Example3/E3_Files.jl")
include("Example3/E3_Results.jl")
include("Example3/sbc.jl")

#include("Sanity.jl")
include("Figures/tau_test_copper.jl")
include("Figures/alpha_test_copper.jl")
include("Figures/aa_nonlip.jl")
include("Figures/aa_test_copper.jl")
include("Figures/aa_alpha_copper.jl")
include("Problem_Data/boundary.jl")
include("Problem_Data/build_problemex1.jl")
include("Problem_Data/Exact.jl")
include("Problem_Data/NL_Equationex1.jl")
include("Example2/NL_Equationex2.jl")
include("Example2/build_problem2.jl")
#include("Example3/build_problem3.jl")

function ExNum(alge::Function)
enum=0
if (alge == alg1e1) || (alge == alg2ge1) || 
                    (alge == GF1!) || (alge == GF2!)
   enum=1
end
if (alge == alg1e2) || (alge == alg2ge2) || 
                    (alge == GF1e2!) || (alge == GF2e2!)
   enum=2
end
if (alge == alg1e3) || (alge == alg2ge3) || 
                    (alge == GF1e3!) || (alge == GF2e3!)
   enum=3
end
return enum
end

export ExNum
export alg1, alg2
export build_problem3, alg1e3, alg2ge3, alg_base3
export alg1_test, alg2g
export alg1e1, uefun_ex12d, alpha_test_copper, tau_test_copper 
export aa_test_copper
export aa_alpha_copper
export testP, aa_test, alg1e1, build_problem, FGen, alg2ge1
export alg1e2, alg2ge2
export complexity, GF1!, GF2!, GF1e2!, GF2e2!, GF1e3!, GF2e3!
export fdgrad, aa_nonlip

end
