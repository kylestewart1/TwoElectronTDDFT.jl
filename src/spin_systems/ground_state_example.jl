include("spin_operators.jl")
include("ground_state.jl")
include("structure_factors.jl")
include("entanglement.jl")


n_qubits = 3
paulis = pauli_operators(n_qubits)

J_perp_list = fill(-1.0, n_qubits-1)
#J_perp_list = [-1.0, 10.0]
J_perp = nearest_neighbor(J_perp_list)

J_par_list = fill(-1.0, n_qubits-1)
J_par = nearest_neighbor(J_par_list)

J_zeros = zeros(size(J_par))

T_op = spin_hamiltonian_base(n_qubits, J_perp, J_zeros)
W_op = spin_hamiltonian_base(n_qubits, J_zeros, J_par)
hz = fill(0.0, n_qubits)
V_op = spin_hamiltonian_add_fields(n_qubits, hz)

H_interacting = T_op + W_op + V_op
H_xy = T_op + V_op


E_xy, ψ_xy = solve_schrodinger_static(H_xy)

# XXZ system
energies, eigenvectors = solve_schrodinger_static(H_interacting)

E₀, ψ₀ = ground_state(E_xy, ψ_xy)

U = rotation_matrix(float(π)/2,[0,1,0], n_qubits)

ψ = U*ψ₀


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

corr_xx = correlation_function(ψ₀, paulis,n_qubits,"x")
corr_yy = correlation_function(ψ₀, paulis,n_qubits,"y")
corr_zz = correlation_function(ψ₀, paulis,n_qubits,"z")


p = LinRange(0.0,pi,n_qubits-1)
struc_xx = structure_factor(p,corr_xx)
struc_yy = structure_factor(p,corr_yy)
struc_zz = structure_factor(p,corr_zz)


coeff_matrix = 2*(J_perp - Diagonal(hz))
λ, O = eigen(coeff_matrix)
eps = 1e-12
λ[abs.(λ) .< eps] .= 0.0
O[abs.(O) .< eps] .= 0.0

λ2, O2 = eigen(H_xy) 
O2[abs.(O2) .< eps] .= 0.0


new_Z = new_pauli_Z(paulis, O)

T = zeros((2^n_qubits, 2^n_qubits))

for k in 1:n_qubits
    H = [0.0 0.0; 0.0 λ[k]]
    for i in 1:(k-1)
        H = kron(Id,H)
    end
    for i in (k+1):n_qubits
        H = kron(H,Id)
    end
    T += H
    print(H)
end

T[abs.(T) .< eps] .= 0.0

eigval, eigvec = eigen(T)
ψ₀_new = eigvec[:,1]

Z_exp = expectation_value(ψ₀, paulis["z"])
Z_exp_new = expectation_value(ψ₀_new, new_Z)