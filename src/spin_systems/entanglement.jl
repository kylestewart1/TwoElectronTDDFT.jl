function two_qubit_concurrence(Z_expectations, k, l)
    N = length(Z_expectations)
    ent = 1 / (N-2)
    for j in [k,l]
        term = sum(Z_expectations) - Z_expectations[j]
        term += (N-3) * Z_expectations[j]
        ent *= sqrt(abs(term))
    end
    return ent
end

function density_matrix(Ψ::Vector{ComplexF64})
    ρ = Ψ * Ψ'
    return ρ
end

function density_matrix(Ψ_t::Vector{Vector{ComplexF64}})
    ρ_t = []
    for i in 1:length(Ψ_t)
        push!(ρ_t, density_matrix(Ψ_t[i]))
    end
    return ρ_t
end

function partial_trace(M::Matrix, which_qubit::Int)
    n_qubits = Int(log(2,size(M)[1]))
    zero, one = [1; 0], [0; 1]
    for i in 1:(which_qubit-1)
        zero = kron(Id,zero)
        one = kron(Id,one)
    end
    for i in (which_qubit+1):n_qubits
        zero = kron(zero,Id)
        one = kron(one,Id)
    end
    partial_tr = zero' * M * zero + one' * M * one
    return partial_tr
end

function partial_trace(M::Matrix, which_qubits::Vector{Int})
    n_qubits = Int(log(2,size(M)[1]))
    which_qubits = sort(which_qubits)
    pt = partial_trace(M, which_qubits[1])
    which_qubits .-= 1
    for i in which_qubits[2:end]
        pt = partial_trace(pt, i)
        which_qubits .-= 1
    end
    return pt
end

function partial_trace(M_t::Vector, which_qubits::Vector{Int})
    pt_t = []
    for i in 1:length(M_t)
        push!(pt_t, partial_trace(M_t[i], which_qubits))
    end
    return pt_t
end


function entropy(ρ::Matrix)
    eigvals = eigen(ρ).values
    S = 0.0
    for eig in eigvals
        if abs(eig) < 1e-10
            eig = 0.0
        end
        if eig != 0
            S -= real(eig*log(eig))
        end
    end
    if abs(S) < 1e-10
        S = 0.0
    end
    return S
end

function entropy(ρ_t::Vector)
    S_t = []
    for i in 1:length(ρ_t)
        push!(S_t, entropy(ρ_t[i]))
    end
    return S_t
end