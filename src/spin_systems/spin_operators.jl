using LinearAlgebra

const Id = [1 0; 0 1]
const X = [0 1; 1 0]
const Y = [0 -im; im 0]
const Z = [1 0; 0 -1]

function single_qubit_gate(n_qubits, which_qubit, gate)
    for j in 1:(which_qubit-1)
        gate = kron(Id, gate)
    end
    for k in (which_qubit+1):n_qubits
        gate = kron(gate, Id)
    end
    return gate    
end


function two_qubit_gate(n_qubits, which_qubits, gates)
    first, second = min(which_qubits...), max(which_qubits...)
    if first == second
        println("Select two different qubits.")
        return nothing
    end
    gate = gates[1]
    for j in 1:(first-1)
        gate = kron(Id, gate)
    end
    for k in (first+1):(second-1)
        gate = kron(gate, Id)
    end
    gate = kron(gate, gates[2])
    for l in (second+1):n_qubits
        gate = kron(gate, Id)
    end
    return gate
end


function current_operator(n_qubits, which_qubits, J_perp)
    first, second = min(which_qubits...), max(which_qubits...)
    if first == second
        println("Select two different qubits.")
        return nothing
    end

    current = two_qubit_gate(n_qubits, (first, second), (X,Y)) - two_qubit_gate(n_qubits, (first, second), (Y,X))
    current *= -2*J_perp
    return current
end


function two_spin_interaction(n_qubits, J, Pauli)
    N_dim = 2^n_qubits 
    H = zeros(N_dim,N_dim)
    for i in 1:(n_qubits-1)
        term = two_qubit_gate(n_qubits, (i,i+1), (Pauli,Pauli))
        term *= J[i]
        H += term
    end
    return H
end

function external_field(n_qubits, h, Pauli=Z)
    N_dim = 2^n_qubits
    H_ext = zeros(N_dim,N_dim)
    for i in 1:n_qubits
        ext_term = single_qubit_gate(n_qubits, i, Pauli)
        ext_term *= h[i]
        H_ext += ext_term
    end
    return H_ext
end


function spin_hamiltonian_base(n_qubits, J_perp, J_par)
    H = two_spin_interaction(n_qubits, J_perp, X) + two_spin_interaction(n_qubits, J_perp, Y)
    H += two_spin_interaction(n_qubits, J_par, Z)
    return H
end

function spin_hamiltonian_add_fields(n_qubits, hz, hx=nothing, hy=nothing)
    H_ext = external_field(n_qubits, hz)
    if hx !== nothing
        H_ext += external_field(n_qubits, hx, X)
    end
    if hy !== nothing
        H_ext += external_field(n_qubits, hy, Y)
    end
    return H_ext
end


function sum_all_but_one(J, gates, excluded_index)
    N = length(gates)
    tot = zeros(size(gates[1]))
    for i in 1:(excluded_index-1)
        tot += J[i]*gates[i]
    end
    for i in (excluded_index+1):N
        tot += J[i] * gates[i]
    end
    return tot
end

function pauli_operators(n_qubits)
    paulis = Dict("x" => [], "y" => [], "z" => [])
    for i in 1:n_qubits
        push!(paulis["x"], single_qubit_gate(n_qubits, i, X))
        push!(paulis["y"], single_qubit_gate(n_qubits, i, Y))
        push!(paulis["z"], single_qubit_gate(n_qubits, i, Z))
    end
    return paulis
end

function current_operators(n_qubits, paulis, J_perp)
    currents = []
    for i in 1:(n_qubits-1)
        push!(currents, -2*J_perp[i] * (paulis["x"][i]*paulis["y"][i+1] 
                                     - paulis["y"][i]*paulis["x"][i+1]))
    end
    return currents
end

function stress_tensor_operator(n_qubits, paulis, J_perp)
    stress_tensor = []
    for i in 1:(n_qubits-1)
        push!(stress_tensor, 4 * J_perp[i] * (
                                       paulis["z"][i]*sum_all_but_one(J_perp,paulis["x"],i)*paulis["x"][i+1]
                                     + paulis["z"][i]*sum_all_but_one(J_perp,paulis["y"],i)*paulis["y"][i+1]
                                     - paulis["x"][i]*sum_all_but_one(J_perp,paulis["x"],i+1)*paulis["z"][i+1]
                                     - paulis["y"][i]*sum_all_but_one(J_perp,paulis["y"],i+1)*paulis["z"][i+1]))
    end
    return stress_tensor
end

function internal_force_density_operator(n_qubits, paulis, J_par)
    internal_force_density = []
    for i in 1:(n_qubits-1)
        push!(internal_force_density, 4 * J_par[i] * (
                                       paulis["y"][i]*sum_all_but_one(J_perp,paulis["z"],i+1)*paulis["y"][i+1]
                                     - paulis["x"][i]*sum_all_but_one(J_perp,paulis["z"],i+1)*paulis["x"][i+1]
                                     + paulis["y"][i]*sum_all_but_one(J_perp,paulis["z"],i)*paulis["y"][i+1]
                                     - paulis["x"][i]*sum_all_but_one(J_perp,paulis["z"],i)*paulis["x"][i+1]))
    end
    return internal_force_density
end

function kinetic_energy_operator(n_qubits, paulis, J_perp)
    kinetic_energy = []
    for i in 1:(n_qubits-1)
        push!(kinetic_energy, J_perp[i] * (paulis["x"][i]*paulis["x"][i+1] + paulis["y"][i]*paulis["y"][i+1]))
    end
    return kinetic_energy
end


function expectation_values(Ψ, operators)
    N = length(operators)
    expectation = []
    for i in 1:N 
        push!(expectation, real(Ψ' * operators[i] * Ψ))
    end
    return expectation
end

function time_dependent_expectation_values(Ψ_t, operators)
    n_time_steps = length(Ψ_t)
    expectation_t = []
    for i in 1:n_time_steps
        expectation = expectation_values(Ψ_t[i], operators)
        push!(expectation_t, expectation)
    end
    return expectation_t
end
