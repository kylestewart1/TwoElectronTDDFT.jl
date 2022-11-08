function solve_schrodinger_static(H)
    eig = eigen(H)
    return eig.values, eig.vectors
end

function ground_state(energies, eigenstates)
    E₀, ψ₀ = minimum(energies), eigenstates[:,argmin(energies)]
    return E₀, ψ₀
end


function heisenberg_energy(Z_exp, J_par, n_qubits)
    E_H = 0.0
    for i in 1:(n_qubits-1)
        E_H += J_par[i,i+1] * Z_exp[i] * Z_exp[i+1]
    end
    return E_H
end

function external_field_energy(Z_exp, h)
    return -sum(h .* Z_exp)
end

function XY_kinetic_energy(X_exp, Y_exp, J_perp)
    T_s = 0.0
    for i in 1:(n_qubits-1)
        T_s += J_perp[i,i+1] * (X_exp[i]X_exp[i+1] + Y_exp[i]Y_exp[i+1])
    end
    return T_s
end

function xc_energy(T, T_s, W, E_H)
    return T - T_s + W - E_H
end


function heisenberg_fields(Z_exp, J_par, n_qubits)
    h_H = zeros(size(Z_exp))
    for i in 1:(n_qubits)
        h = 0.0
        if i > 1
            h += J_par[i-1,i] * Z_exp[i-1]
        end
        if i < n_qubits
            h += J_par[i,i+1] * Z_exp[i+1]
        end
        h_H[i] = h
    end
    return -h_H
end