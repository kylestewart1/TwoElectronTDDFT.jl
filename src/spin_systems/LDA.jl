include("spin_operators.jl")
include("ground_state.jl")


function homogeneous_state(uniform_Z, n_qubits)
    Ψ = zeros(2^n_qubits)
    Ψ[1] = 1.0
    θ = acos(uniform_Z)
    U = rotation_matrix(θ, [0,1,0], n_qubits)
    return U*Ψ
end

function energy_functional(n_qubits::Int, J_perp::Matrix, J_par::Matrix, uniform_Z::Float64)
    H = spin_hamiltonian_base(n_qubits, J_perp, J_par)
    Ψ = homogeneous_state(uniform_Z, n_qubits)
    E = expectation_value(Ψ, H)
    return E
end

function energy_functional(n_qubits::Int, J_perp::Matrix, J_par::Matrix, uniform_Z::LinRange)
    E = zeros(size(uniform_Z))
    for i in eachindex(uniform_Z)
        E[i] = energy_functional(n_qubits, J_perp, J_par, uniform_Z[i])
    end
    return E
end

function energy_functional(n_qubits::Int, J_perp::Matrix, scaling::Vector, uniform_Z)
    E = []
    for i in eachindex(scaling)
        J_par = scaling[i] .* J_perp
        push!(E, energy_functional(n_qubits,J_perp,J_par,uniform_Z))
    end
    return E
end


function energy_functional(n_qubits::Vector, uniform_Z)
    E = []
    scaling = [-10.0, -1.0, 0.0, 1.0, 10.0, 100.0]
    for i in eachindex(n_qubits)
        J_perp_list = fill(-1.0, n_qubits[i]-1)
        J_perp = nearest_neighbor(J_perp_list)
        push!(E, energy_functional(n_qubits[i], J_perp, scaling, uniform_Z))
    end
    return E
end



function LDA(σᶻ::Float64, J_perp_list::Vector, J_par_list::Vector)
    γ = sum(J_perp_list)
    λ = sum(J_par_list)
    return γ + (λ - γ)*σᶻ^2
end

function  LDA(σᶻ::LinRange, J_perp_list::Vector, J_par_list::Vector)
    γ = sum(J_perp_list)
    λ = sum(J_par_list)
    return γ .+ (λ - γ)*σᶻ.^2
end

function LDA(Z_list::Vector, J_perp_list::Vector, J_par_list::Vector)
    γ = sum(J_perp_list)
    λ = sum(J_par_list)
    n_qubits = length(J_perp_list)+1
    E = γ
    """
    for i in eachindex(Z_list)
        E += ((λ - γ) / n_qubits) * Z_list[i]^2
    end
    """
    s = sum(Z_list) / n_qubits
    E += (λ - γ) * s^2
    return E
end

function xy_kinetic_plus_correlation(Z_list::Vector, J_perp_list::Vector, J_par_list::Vector, h_ext::Vector)
    E_tot = LDA(Z_list, J_perp_list, J_par_list)
    E_ext = external_field_energy(Z_list, h_ext)
    E_H = heisenberg_energy(Z_list, J_par_list)
    Ts_plus_Ec = E_tot - E_ext - E_H
    return Ts_plus_Ec
end

