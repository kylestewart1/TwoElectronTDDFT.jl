include("propagation.jl")

n_qubits = 3
Δt = 0.01
n_time_steps = 500
t = 0:Δt:n_time_steps*Δt
hz = [cos.(t) cos.(t .- 2*pi/3) cos.(t .- 4*pi/3)]
Ψ₀ = ComplexF64[0.0; 1.0; 1.0; 0.0; 1.0; 0.0; 0.0; 0.0]
Ψ₀ /= norm(Ψ₀)
J_perp = [1.0 1.0 1.0]
J_par = [1.0 1.0 1.0]

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



# inversion

include("inversion.jl")

paulis = pauli_operators(n_qubits)
current_ops = current_operators(n_qubits, paulis, J_perp)
currents = time_dependent_expectation_values(Ψ_t, current_ops)

current_derivs = get_current_derivs(currents,Δt)

J_par_prime = copy(J_par) #[0.0 0.0 0.0]
J_perp_prime = copy(J_perp)

Ψ₀_prime = copy(Ψ₀)
Z₀_prime = expectation_values(Ψ₀_prime, Z_operators)
H_base_prime = spin_hamiltonian_base(n_qubits, J_perp_prime, J_par_prime)

global_field = copy(hz[:,1])



Ψ_t_prime, h_t_prime = inversion(Ψ₀_prime, current_derivs, n_qubits, J_perp_prime, J_par_prime, Δt, global_field)
Z_t_prime = time_dependent_expectation_values(Ψ_t_prime, Z_operators)
Z_1_prime = Float64[]
Z_2_prime = Float64[]
Z_3_prime = Float64[]
for i in 1:length(Z_t_prime)
    push!(Z_1_prime, Z_t_prime[i][1])
    push!(Z_2_prime, Z_t_prime[i][2])
    push!(Z_3_prime, Z_t_prime[i][3])
end

z1_prime_plot = plot(t, Z_1_prime, label="qubit 1", xlabel="Time", ylabel="Z expectation value",
    title="Z expectation values in auxiliary system")
plot!(t, Z_2_prime, label="qubit 2", )
plot!(t, Z_3_prime, label="qubit 3")

hz_prime = reduce(hcat,h_t_prime)'

h_prime_plot = plot(t[1:(end-1)], hz_prime[:,1], label="qubit 1", xlabel = "Time", ylabel = "External field", 
    title="External fields in Z direction")
plot!(t[1:(end-1)], hz_prime[:,2], label="qubit 2", )
plot!(t[1:(end-1)], hz_prime[:,3], label="qubit 3")
aux_plt = plt = plot(h_prime_plot, z1_prime_plot, layout=(2,1))

savefig(aux_plt, "Z_and_fields_auxiliary_system.png")
