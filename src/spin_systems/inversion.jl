function two_point_difference(f, f_next, Δt)
    return (f_next .- f) ./ Δt
end

function three_point_midpoint(f::Array, Δt)
    return (f[3]-f[1])/(2*Δt)
end

function five_point_left_endpoint(f::Array, Δt)
    return (-25*f[1] + 48*f[2]-36*f[3]+16*f[4]-3*f[5])/(12*Δt)
end

five_point_right_endpoint(f::Array, Δt) = -five_point_left_endpoint(f[end:-1:1],Δt)

function five_point_midpoint(f::Array,Δt)
    return (f[1] - 8*f[2] + 8*f[4] - f[5])/(12*Δt)
end



function get_current_derivs(currents, Δt)
    derivs = []
    n_time_steps = length(currents)
    for n in 1:(n_time_steps-1)
        if n <= 2
            deriv = five_point_left_endpoint(currents[n:n+5], Δt)
        elseif n >= (n_time_steps-2)
            deriv = five_point_right_endpoint(currents[n-4:n], Δt)
        else
            deriv = five_point_midpoint(currents[n-2:n+2], Δt)
        end
        push!(derivs, deriv)
    end
    return derivs
end


function solve_for_fields(current_derivs, stress_tensor, internal_force_density, kinetic_energy, fixed_field, n_qubits, fixed_field_index=1)
    h = zeros(n_qubits)
    h[fixed_field_index] = copy(fixed_field)
    δh = zeros(n_qubits-1)
    for i in 1:(n_qubits-1)
        δh[i] = (current_derivs[i,i+1] - stress_tensor[i,i+1] - internal_force_density[i,i+1]) / (4*kinetic_energy[i,i+1])
        if abs(δh[i]) < 1e-5
            δh[i]=0
        end
    end
    for i in 1:(n_qubits-1)
        h[i+1] = h[i] - δh[i]
    end
    return h
end

function update_hamiltonian(n_qubits, H_base, h)
    H = H_base + spin_hamiltonian_add_fields(n_qubits, h)
    return H
end


function advance_wf(Ψ, H, Δt)
    Ψ_next = Ψ - (H * Ψ) .* Δt * im
    return Ψ_next
end

function inversion_time_step(Ψ, current_derivs, stress_tensor_op, internal_force_op, kinetic_op, fixed_field, n_qubits, Δt, H_base)
    stress_tensor = expectation_value(Ψ, stress_tensor_op)
    internal_force_density = expectation_value(Ψ, internal_force_op)
    kinetic_energy = expectation_value(Ψ, kinetic_op)
    h = solve_for_fields(current_derivs, stress_tensor, internal_force_density, kinetic_energy, fixed_field, n_qubits)
    H = update_hamiltonian(n_qubits, H_base, h)
    """
    Ψ_mid = advance_wf(Ψ, H, Δt/2)
    stress_tensor_mid = expectation_value(Ψ_mid, stress_tensor_op)
    internal_force_density_mid = expectation_value(Ψ_mid, internal_force_op)
    kinetic_energy_mid = expectation_value(Ψ_mid, kinetic_op)
    h_mid = solve_for_fields(current_derivs, stress_tensor_mid, internal_force_density_mid, kinetic_energy_mid, fixed_field, n_qubits)
    H_mid = update_hamiltonian(n_qubits, H_base, h_mid)
    Ψ_next = crank_nicholson(Ψ, H_mid, Δt)
    """
    Ψ_next = advance_wf(Ψ, H, Δt)
    Ψ_next /= norm(Ψ_next)
    return h, Ψ_next
end
   

function inversion(Ψ₀, current_derivs, n_qubits, J_perp, J_par, Δt, fixed_field_t)
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
        h, Ψ = inversion_time_step(Ψ, current_derivs[i], stress_tensor_op, internal_force_op, kinetic_op, fixed_field_t[i], n_qubits, Δt, H_base)
        push!(Ψ_t, Ψ)
        push!(h_t, h) 
    end
    return Ψ_t, h_t
end
