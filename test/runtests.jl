using TwoElectronTDDFT
using Test

@testset "Operators" begin
    x = -1.0:1:1.0
    V = Potential(soft_coulomb(x))
    display(V.data)
end
