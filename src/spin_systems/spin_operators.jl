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
    gate_1 = single_qubit_gate(n_qubits, which_qubits[1], gates[1]) 
    gate_2 = single_qubit_gate(n_qubits, which_qubits[2], gates[2]) 
    gate = gate_1*gate_2
    return gate
end


function current_operator(n_qubits, which_qubits, J_perp)
    current = two_qubit_gate(n_qubits, which_qubits, (X,Y)) - two_qubit_gate(n_qubits, which_qubits, (Y,X))
    current *= -2*J_perp
    return current
end


function two_spin_interaction(n_qubits, J, Pauli)
    N_dim = 2^n_qubits 
    H = zeros(N_dim,N_dim)
    for i in 1:n_qubits, j in (i+1):n_qubits
        term = two_qubit_gate(n_qubits, (i,j), (Pauli,Pauli))
        term *= J[i,j]
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
    return -H_ext
end


function sum_all_but_one(J, gates, excluded_index)
    N = length(gates)
    tot = zeros(size(gates[1]))
    for m in 1:(excluded_index-1)
        tot += J[m,excluded_index]*gates[m]
    end
    for m in (excluded_index+1):N
        tot += J[m,excluded_index] * gates[m]
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
    current_ops = zeros(ComplexF64,(n_qubits,n_qubits,2^n_qubits,2^n_qubits))
    for i in 1:n_qubits, j in (i+1):n_qubits
        current_ops[i,j,:,:] = current_operator(n_qubits,(i,j),J_perp[i,j])
    end
    for i in 1:n_qubits, j in 1:i
        current_ops[i,j,:,:] = -current_ops[j,i,:,:]
    end
    return current_ops
end

function stress_tensor_operator(n_qubits, paulis, J_perp)
    stress_tensor = zeros((n_qubits,n_qubits,2^n_qubits,2^n_qubits))
    for k in 1:n_qubits, i in (k+1):n_qubits
        zx = -paulis["z"][k]*sum_all_but_one(J_perp,paulis["x"],k)*paulis["x"][i]
        zy = -paulis["z"][k]*sum_all_but_one(J_perp,paulis["y"],k)*paulis["y"][i]
        xz = paulis["x"][k]*sum_all_but_one(J_perp,paulis["x"],i)*paulis["z"][i]
        yz = paulis["y"][k]*sum_all_but_one(J_perp,paulis["y"],i)*paulis["z"][i]
        stress_tensor[k,i,:,:] = 4*J_perp[k,i]*(zx + zy + xz + yz)
    end
    for k in 1:n_qubits, i in 1:k
        stress_tensor[k,i,:,:] = -stress_tensor[i,k,:,:]
    end
    return stress_tensor
end

function internal_force_density_operator(n_qubits, paulis, J_par)
    internal_force_density = zeros((n_qubits,n_qubits,2^n_qubits,2^n_qubits))
    for k in 1:n_qubits, i in (k+1):n_qubits
        yz1 = paulis["y"][k]*sum_all_but_one(J_par,paulis["z"],i)*paulis["y"][i]
        xz1 = - paulis["x"][k]*sum_all_but_one(J_par,paulis["z"],i)*paulis["x"][i]
        yz2 = paulis["y"][k]*sum_all_but_one(J_par,paulis["z"],k)*paulis["y"][i]
        xz2 = -paulis["x"][k]*sum_all_but_one(J_par,paulis["z"],k)*paulis["x"][i]
        internal_force_density[k,i,:,:] = 4*J_par[k,i]*(yz1 + xz1 + yz2 + xz2)
    end
    for k in 1:n_qubits, i in 1:k
        internal_force_density[k,i,:,:] = -internal_force_density[i,k,:,:]
    end
    return internal_force_density
end

function kinetic_energy_operator(n_qubits, paulis, J_perp)
    kinetic_energy = zeros(n_qubits,n_qubits,2^n_qubits,2^n_qubits)
    for i in 1:n_qubits, j in i:n_qubits
        kinetic_energy[i,j,:,:] = J_perp[i,j] * (paulis["x"][i]*paulis["x"][j] + paulis["y"][i]*paulis["y"][j])
    end
    for i in 1:n_qubits, j in 1:(i-1)
        kinetic_energy[i,j,:,:] = kinetic_energy[j,i,:,:]
    end
    return kinetic_energy
end


function expectation_value(Ψ, operator::Matrix)
    return real(Ψ' * operator * Ψ)
end

function expectation_value(Ψ, operators::Vector)
    expectation = []
    for n in 1:length(operators)
        push!(expectation, expectation_value(Ψ, operators[n]))
    end
    return expectation
end


function expectation_value(Ψ, operators::Array)
    op_shape = size(operators)[1:2]
    expectation = zeros(op_shape)
    operators = complex(operators)
    for i in 1:op_shape[1], j in 1:op_shape[2]
        expectation[i,j] = expectation_value(Ψ, operators[i,j,:,:])
    end
    return expectation
end

function time_dependent_expectation_values(Ψ_t, operators)
    n_time_steps = length(Ψ_t)
    expectation_t = []
    for i in 1:n_time_steps
        expectation = expectation_value(Ψ_t[i], operators)
        push!(expectation_t, expectation)
    end
    return expectation_t 
end

function time_dependent_expectation_values(Ψ_t, operator::Matrix)
    n_time_steps = length(Ψ_t)
    expectation_t = []
    for i in 1:n_time_steps
        expectation = expectation_value(Ψ_t[i], operator)
        push!(expectation_t, expectation)
    end
    return expectation_t 
end


function time_dependent_expectation_values(Ψ_t, operators::Array{Number, 4})
    n_time_steps = length(Ψ_t)
    op_shape = size(operators)[1:2]
    expectation_t = []
    for n in 1:n_time_steps
        expectation = expectation_value(Ψ_t[n], operators)
        push!(expectation_t, expectation)
    end
    return expectation_t 
end


function hadamard_operator(n_qubits)
    H1 = (1/sqrt(2)) * [1 1; 1 -1]
    H = copy(H1)
    for i in 1:(n_qubits-1)
        H =  kron(H, H1)
    end
    return H
end

function nearest_neighbor(J_list::Vector)
    n_qubits = size(J_list)[1]+1
    J = zeros((n_qubits,n_qubits))
    for i in 1:(n_qubits-1)
        J[i,i+1] = J_list[i]
        J[i+1,i] = J_list[i]
    end
    return J
end

function rotation_matrix(θ::Float64, n::Vector)
    n ./= norm(n)
    R = cos(θ/2)*I - im*sin(θ/2)*(n[1]*X + n[2]*Y + n[3]*Z)
    return R
end

function rotation_matrix(θ::Float64, n::Vector, n_qubits::Int)
    R_1 = rotation_matrix(θ, n)
    R_N = R_1
    for i in 1:(n_qubits-1)
        R_N = kron(R_1, R_N)
    end
    return R_N
end

import Base.round
function round(A::Array{ComplexF64}, ϵ=1e-8)
    B = copy(A)
    for i in eachindex(B)
        x, y = real(B[i]), imag(B[i])
        if abs(x) < ϵ
            x = 0.0
        end
        if abs(y) < ϵ
            y = 0.0
        end
        B[i] = x + y*im
    end
    return B
end

function pauli_change_basis(paulis, M::AbstractMatrix)
    new_paulis = Dict("x" => [], "y" => [], "z" => [])
    for α in eachindex(paulis)
        for i in eachindex(paulis[α])
            push!(new_paulis[α], sum(M[i,:] .* paulis[α]))
        end
    end
    return new_paulis
end
"""
function new_pauli_Z(paulis, M::AbstractMatrix)
    new_Zs = []
    n_qubits = length(paulis["z"])
    for i in eachindex(paulis["z"])
        Z_new = zeros(size(paulis["z"][1]))
        for k in 1:n_qubits, l in 1:n_qubits
            Z_new += (-1/2) * M[i,k]*conj(M[i,l])*(paulis["x"][k]*paulis["x"][l] 
                    + paulis["y"][k]*paulis["y"][l])
        end
        push!(new_Zs, Z_new)
    end
    return new_Zs
end
"""

function new_pauli_Z(paulis, U::AbstractMatrix)
    new_Zs = []
    n_qubits = length(paulis["z"])
    for i in eachindex(paulis["z"])
        Z_new = zeros(size(paulis["z"][1]))
        for l in 1:n_qubits
            for k in 1:(l-1)
                term = (-1/2) * U[i,k]*conj(U[i,l])*(paulis["x"][k] + im*paulis["y"][k])
                for m in (k+1):(l-1)
                    term *= paulis["z"][m]
                end
                term *= (paulis["x"][l] - im*paulis["y"][l])
                Z_new += term
            end
            for k in (l+1):n_qubits
                term = (-1/2) * U[i,k]*conj(U[i,l])*(paulis["x"][l] - im*paulis["y"][l])
                for m in (l+1):(k-1)
                    term *= paulis["z"][m]
                end
                term *= (paulis["x"][k] + im*paulis["y"][k])
                Z_new += term
            end
        end
        for k in 1:n_qubits
            Z_new += U[i,k]*conj(U[i,k])*paulis["z"][k]
        end

        push!(new_Zs, Z_new)
    end
    return new_Zs
end


function swap_col!(M::Matrix, i, j)
    for k in axes(M)[1]
        idata = M[k,i]
        M[k,i] = M[k,j]
        M[k,j] = idata
    end
end