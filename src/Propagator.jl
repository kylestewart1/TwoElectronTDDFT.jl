struct CrankNicholson <: AbstractPropagator
    dt::Float64
    H_mid::Hamiltonian
end

function (U::CrankNicholson)(ψ::Ket)
    Δ = U.H_mid * (im/2.0) * U.dt
    lhs = I .+ Δ
    rhs = (I .- Δ) * ψ    
    ψ_new = lhs / rhs # solves the linear equation lhs*X = rhs
    update!(ψ, ψ_new)
end


function update!(U::CrankNicholson, H_mid::Hamiltonian)
    for i in eachindex(H_mid)
        U.H_mid[i] = H_mid[i]
    end
end