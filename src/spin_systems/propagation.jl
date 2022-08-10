include("spin_operators.jl")
include("entanglement.jl")

function crank_nicholson(Ψ, H_mid, Δt)
    lhs = I + H_mid * (Δt * im /2)
    rhs = (I - H_mid * (Δt * im /2)) * Ψ
    Ψ_new = lhs \ rhs
    return Ψ_new
end


function time_evolution(Ψ₀, Δt, n_time_steps, H_base, hz_ext)
    Ψ_t = [Ψ₀]
    for i in 1:(n_time_steps)
        hz_mid = (hz_ext[i,:] + hz_ext[i+1,:])/2
        H_mid = H_base + spin_hamiltonian_add_fields(n_qubits, hz_mid)
        Ψ_new = crank_nicholson(Ψ_t[i], H_mid, Δt)
        push!(Ψ_t, Ψ_new)
    end
    return Ψ_t
end