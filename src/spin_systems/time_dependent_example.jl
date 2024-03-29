include("propagation.jl")
include("ground_state.jl")

n_qubits = 3
Δt = 0.0001
n_time_steps = 50000
t = 0:Δt:n_time_steps*Δt
hz = [cos.(t) cos.(t .- 2*pi/3) cos.(t .- 4*pi/3)]
#hz = zeros(length(t),n_qubits)

# 3 qubit example
Ψ₀ = ComplexF64[0.0; 1.0; 1.0; 0.0; 1.0; 0.0; 0.0; 0.0]
Ψ₀ /= norm(Ψ₀)
J_perp = [0.0 -1.0 0.0; -1.0 0.0 -1.0; 0.0 -1.0 0.0]
J_par = [0.0 -0.01 0.0; -0.01 0.0 -0.01; 0.0 -0.01 0.0]

# 2 qubit example
"""
Ψ₀ = ComplexF64[0.0; 1.0; 1.0; 0.0]
Ψ₀ /= norm(Ψ₀)
J_perp = [0.0 1.0; 1.0 0.0]
J_par = [0.0 1.0; 1.0 0.0]
"""

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
plot!(t, hz[:,2], label="qubit 2")
plot!(t, hz[:,3], label="qubit 3")
    
plt = plot(h_plot, z1_plot, layout=(2,1))
savefig(plt, "original_system.png")


# inversion

include("inversion.jl")

paulis = pauli_operators(n_qubits)
current_ops = current_operators(n_qubits, paulis, J_perp)
currents = time_dependent_expectation_values(Ψ_t, current_ops)

current_1_2 = Float64[]
current_2_3 = Float64[]
for i in 1:length(currents)
    push!(current_1_2, currents[i][1,2])
    push!(current_2_3, currents[i][2,3])
end
currents_plot = plot(t, current_1_2, label="j_1_2", xlabel = "Time", ylabel="Current expectation value")
plot!(t, current_2_3, label="j_2_3")

current_derivs = get_current_derivs(currents,Δt)

J_perp_prime = [0.0 -1.0 0.0; -1.0 0.0 -1.0; 0.0 -1.0 0.0]
J_par_prime = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]


Ψ₀_prime = copy(Ψ₀)
Z₀_prime = expectation_value(Ψ₀_prime, Z_operators)
H_base_prime = spin_hamiltonian_base(n_qubits, J_perp_prime, J_par_prime)

#global_field = copy(hz[:,1])
global_field = zeros(n_time_steps+1)


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
plot!(t, Z_2_prime, label="qubit 2")
plot!(t, Z_3_prime, label="qubit 3")

hz_prime = reduce(hcat,h_t_prime)'

h_prime_plot = plot(t[1:(end-1)], hz_prime[:,1], label="qubit 1", xlabel = "Time", ylabel = "External field", 
    title="External fields in Z direction")
plot!(t[1:(end-1)], hz_prime[:,2], label="qubit 2")
plot!(t[1:(end-1)], hz_prime[:,3], label="qubit 3")
aux_plt = plt = plot(h_prime_plot, z1_prime_plot, layout=(2,1))
savefig(aux_plt, "auxiliary_system.png")


Z_1_diff = Z_1_prime - Z_1
Z_2_diff = Z_2_prime - Z_2
Z_3_diff = Z_3_prime - Z_3
Z_diff_plt = plot(t, Z_1_diff, label="qubit 1", xlabel="Time", ylabel="Z expectation value",
title="Difference in Z")
plot!(t, Z_2_diff, label="qubit 2")
plot!(t, Z_3_diff, label="qubit 3")
original_vs_aux_plt = plot(z1_plot, z1_prime_plot, Z_diff_plt, layout=(3,1))
savefig("z_difference.png")



ρ_t = density_matrix(Ψ_t)

ρ_1_t = partial_trace(ρ_t, [2,3])
ρ_2_t = partial_trace(ρ_t, [1,3])
ρ_3_t = partial_trace(ρ_t, [1,2])

S_1_t = entropy(ρ_1_t)
S_2_t = entropy(ρ_2_t)
S_3_t = entropy(ρ_3_t)

entropy_plt = plot(t, S_1_t, label="qubit 1", xlabel="Time", ylabel="Entanglement entropy",
    title="Entropy of entanglement")
plot!(t, S_2_t, label="qubit 2")
plot!(t, S_3_t, label="qubit 3")


ρ_t_prime = density_matrix(Ψ_t_prime)

ρ_1_t_prime = partial_trace(ρ_t_prime, [2,3])
ρ_2_t_prime = partial_trace(ρ_t_prime, [1,3])
ρ_3_t_prime = partial_trace(ρ_t_prime, [1,2])

S_1_t_prime = entropy(ρ_1_t_prime)
S_2_t_prime = entropy(ρ_2_t_prime)
S_3_t_prime = entropy(ρ_3_t_prime)

entropy_plt_prime = plot(t, S_1_t_prime, label="qubit 1", xlabel="Time", ylabel="Entanglement entropy",
    title="Auxiliary system")
plot!(t, S_2_t_prime, label="qubit 2")
plot!(t, S_3_t_prime, label="qubit 3")

diff_1 = S_1_t_prime - S_1_t
diff_2 = S_2_t_prime - S_2_t
diff_3 = S_3_t_prime - S_3_t

diff_plt = plot(t, diff_1, label="qubit 1", xlabel="Time", ylabel="Entanglement entropy",
title="Difference")
plot!(t, diff_2, label="qubit 2")
plot!(t, diff_3, label="qubit 3")


entanglement_plt = plot(entropy_plt, entropy_plt_prime, diff_plt, layout=(3,1))
savefig(entropy_plt, "entanglement_entropy.png")
