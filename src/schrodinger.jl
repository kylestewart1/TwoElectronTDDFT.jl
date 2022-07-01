include("spin_operators.jl")
include("density.jl")
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

n_qubits = 3
Δt = 0.01
n_time_steps = 601
t = 0:Δt:n_time_steps*Δt
hz = [cos.(t) cos.(t .- 2*pi/3) cos.(t .- 4*pi/3)]
Ψ₀ = ComplexF64[0.0; 1.0; 1.0; 0.0; 1.0; 0.0; 0.0; 0.0]
Ψ₀ /= norm(Ψ₀)
J_perp = [1 1 1]
J_par = [1 1 1]

H_base = spin_hamiltonian_base(n_qubits, J_perp, J_par)

Ψ_t = time_evolution(Ψ₀, Δt, n_time_steps, H_base, hz)

Z_operators = [single_qubit_gate(n_qubits, i, Z) for i in 1:n_qubits]

Z_t = time_dependent_expectation_values(Ψ_t, Z_operators)


using Plots
Z_1 = Float64[]
Z_2 = Float64[]
Z_3 = Float64[]
for i in eachindex(t)
    push!(Z_1, Z_t[i][1])
    push!(Z_2, Z_t[i][2])
    push!(Z_3, Z_t[i][3])
end
z1_plot = plot(t, Z_1, label="qubit 1", xlabel="Time", ylabel="Z expectation value",
    title="Z expectation values")
plot!(t, Z_2, label="qubit 2", )
plot!(t, Z_3, label="qubit 3")

h_plot = plot(t, hz[:,1], label="qubit 1", xlabel = "Time", ylabel = "External field", 
    title="External fields in Z direction")
plot!(t, hz[:,2], label="qubit 2", )
plot!(t, hz[:,3], label="qubit 3")
    
plt = plot(h_plot, z1_plot, layout=(2,1))
savefig(plt, "Z_and_fields_plot.png")

entanglement = [Float64[], Float64[], Float64[]]
for i in eachindex(t)
    push!(entanglement[1], two_qubit_entanglement(Z_t[i], 1, 2))
    push!(entanglement[2], two_qubit_entanglement(Z_t[i], 1, 3))
    push!(entanglement[3], two_qubit_entanglement(Z_t[i], 2, 3))
end

ent_plt = plot(t, entanglement[1], label="qubits 1 and 2", xlabel="Time", ylabel="Concurrence", 
title = "Two-qubit concurrence")
plot!(t, entanglement[2], label="qubits 1 and 3")
plot!(t, entanglement[3], label="qubits 2 and 3")
savefig(ent_plt, "concurrence.png")