using LinearAlgebra


function mean_field_energy(spins::Vector, J::Matrix)
    E = 0.0
    for i in eachindex(spins), j in (i+1):size(spins)[1]
        E += J[i,j] * (spins[i]' * spins[j])
    end
    return E
end

mean_field_energy(spins::Vector, J::Vector) = mean_field_energy(spins, diagm(1 => J, -1 => J))
mean_field_energy(spins::Vector, J::Float64) = mean_field_energy(spins, fill(J,size(spins)[1]-1))

function external_field_energy(spins::Vector, B::Vector{<:Real})
    E = 0.0
    for i in eachindex(spins)
        E += spins[i]' * B
    end
    return E
end

function external_field_energy(spins::Vector, B::Vector{<:AbstractVector})
    E = 0.0
    for i in eachindex(spins)
        E += spins[i]' * B[i]
    end
    return E
end

function lsa_heisenberg_afm_1d(spins::Vector, J::Float64, S::Float64)
    E = 0.0
    for i in eachindex(spins)
        E += norm(spins[i])
    end
    E *= J * (2/pi -1) * sqrt(S / (S+1))
    return E
end

function lba_heisenberg_afm_1d(S::Float64, J::Matrix)
    J_unique = UpperTriangular(J)
    N_bonds = count(!=(0.0),J_unique)
    N_spins = size(J)[1]
    E = sum(J_unique)
    E *= S * (2/pi - 1) * N_spins / N_bonds
    return E
end

lba_heisenberg_afm_1d(S::Float64, J::Vector) = lba_heisenberg_afm_1d(S, diagm(1 => J))


function lbsa_heisenberg_afm_1d(spins::Vector, J::Matrix, S::Float64)
    J_unique = UpperTriangular(J)
    N_bonds = count(!=(0.0), J_unique)
    spin_factor = sum(norm.(spins))
    bond_factor = sum(J_unique)
    E = spin_factor * bond_factor * (2/pi - 1) * sqrt(S / (S+1)) / N_bonds
    return E
end


function total_energy(spins, J, B, correlation)
    E = mean_field_energy(spins, J) + external_field_energy(spins, B) + correlation(spins, J)
    return E
end

function spin_constraint_L1(spin_vectors, S)
    diff = norm.(spin_vectors).^2 .- S*(S+1)
    return sum(abs.(diff))
end

function spin_constraint_L2(spin_vectors, S)
    diff = norm.(spin_vectors).^2 .- S*(S+1)
    return sum(diff.^2)
end


function functional_with_constraint(spins, J, B, correlation, constraint, S, λ)
    return total_energy(spins, J, B, correlation) + λ * constraint(spins, S)
end


function normalize_spin(spin::Vector, S::Float64)
    spin .*= sqrt(S*(S+1)) ./ norm.(spin)
    return spin
end