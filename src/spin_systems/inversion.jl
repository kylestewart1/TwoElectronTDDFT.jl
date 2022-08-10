function two_point_difference(f, f_next,Δt)
    return (f_next - f) / Δt
end

function get_current_derivs(currents, Δt)
    current_derivs = []
    n_time_steps = length(currents)
    n_currents = length(currents[1])
    for i in 1:(n_time_steps-1)
        derivs = []
        for j in 1:n_currents
            push!(derivs, two_point_difference(currents[i][j], currents[i+1][j], Δt))
        end
        push!(current_derivs, derivs)
    end
    return current_derivs
end


function solve_for_fields(current_derivs, stress_tensor, internal_force_density, kinetic_energy, global_field)
    δh = []
    for i in 1:length(current_derivs)
        push!(δh, (current_derivs[i] - stress_tensor[i] - internal_force_density[i]) / (4*kinetic_energy[i]))
    end
    h = zeros(length(δh)+1)
    h[1] = copy(global_field)
    for i in 1:(length(h)-1)
        h[i+1] = h[i] - δh[i]
    end
    return h
end

function update_hamiltonian(n_qubits, H_base, h)
    H = H_base + spin_hamiltonian_add_fields(n_qubits, h)
    return H
end


function advance_wf(Ψ, H, Δt)
    Ψ_next = Ψ +  (H * Ψ) .* Δt
    return Ψ_next
end

function inversion_time_step(Ψ, current_derivs, stress_tensor_op, internal_force_op, kinetic_op, global_field, Δt, H_base)
    stress_tensor = expectation_values(Ψ, stress_tensor_op)
    internal_force_density = expectation_values(Ψ, internal_force_op)
    kinetic_energy = expectation_values(Ψ, kinetic_op)
    h = solve_for_fields(current_derivs, stress_tensor, internal_force_density, kinetic_energy, global_field)
    H = update_hamiltonian(n_qubits, H_base, h)
    Ψ_next = advance_wf(Ψ, H, Δt)
    Ψ_next /= norm(Ψ_next)
    return h, Ψ_next
end
        

function inversion(Ψ₀, current_derivs, n_qubits, J_perp, J_par, Δt, global_field_t)
    n_time_steps = length(current_derivs)
    paulis = pauli_operators(n_qubits)
    stress_tensor_op = stress_tensor_operator(n_qubits, paulis, J_perp)
    internal_force_op = internal_force_density_operator(n_qubits, paulis, J_par)
    kinetic_op = kinetic_energy_operator(n_qubits, paulis, J_perp)
    H_base = spin_hamiltonian_base(n_qubits, J_perp, J_par)
    Ψ_t = [Ψ₀]
    h_t = []
    Ψ = Ψ₀
    for i in 1:n_time_steps
       h, Ψ = inversion_time_step(Ψ, current_derivs[i], stress_tensor_op, internal_force_op, kinetic_op, global_field_t[i], Δt, H_base)
       push!(Ψ_t, Ψ)
       push!(h_t, h) 
    end
    return Ψ_t, h_t
end
