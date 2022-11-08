include("spin_operators.jl")
include("ground_state.jl")

n_qubits = 3
paulis = pauli_operators(n_qubits)
J_perp = [0.0 -1.0 0.0; -1.0 0.0 -1.0; 0.0 -1.0 0.0]
J_par = [0.0 -10.0 0.0; -10.0 0.0 -10.0; 0.0 -10.0 0.0]
J_zeros = zeros(size(J_par))

T_op = spin_hamiltonian_base(n_qubits, J_perp, J_zeros)
W_op = spin_hamiltonian_base(n_qubits, J_zeros, J_par)
hz = [10.0, 10.0, 10.0]
V_op = spin_hamiltonian_add_fields(n_qubits, hz)

H_interacting = T_op + W_op + V_op


energies, eigenvectors = solve_schrodinger_static(H_interacting)

E₀, ψ₀ = ground_state(energies, eigenvectors)

X_exp = expectation_value(ψ₀, paulis["x"])
Y_exp = expectation_value(ψ₀, paulis["y"])
Z_exp = expectation_value(ψ₀, paulis["z"])


E_ext = external_field_energy(Z_exp, hz)
T_s = XY_kinetic_energy(X_exp, Y_exp, J_perp)
T = expectation_value(ψ₀, T_op)
E_H = heisenberg_energy(Z_exp, J_par, n_qubits)
W = expectation_value(ψ₀, W_op)
E_xc = xc_energy(T, T_s, W, E_H)



h_H = heisenberg_fields(Z_exp, J_par, n_qubits)