struct CrankNicholson <: Propagator
    dt::Float64
    H_mid::Hamiltonian
end

function (U::CrankNicholson)(ψ::Ket)
    Δ = U.H_mid * (im/2.0) * U.dt
    lhs = I .+ Δ
    rhs = (I .- Δ) * ψ    
    lhs / rhs # solves the linear equation lhs*X = rhs
end