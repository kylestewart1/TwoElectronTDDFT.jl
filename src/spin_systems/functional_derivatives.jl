function mean_field_deriv(spins::Vector, J::Matrix)
    derivs = []
    for i in eachindex(spins)
        deriv = zeros(3)
        for j in eachindex(spins)
            if j != i
                deriv += J[i,j] * spins[j]
            end
        end
        push!(derivs,deriv) 
    end
    return derivs
end

mean_field_deriv(spins::Vector, J::Vector) = mean_field_deriv(spins, diagm(1 => J, -1 => J))
mean_field_deriv(spins::Vector, J::Float64) = mean_field_deriv(spins, fill(J,size(spins)[1]-1))

function lsa_heisenberg_derivs(spins::Vector, J)
    derivs = []
    prefactor = J * (2/pi - 1) * sum(norm.(spins))
    for i in eachindex(spins)
        deriv = prefactor*spins[i]
        push!(derivs, deriv)
    end
    return derivs
end


function spin_constraint_L2_derivs(spins::Vector, S::Float64)
    derivs = []
    for i in eachindex(spins)
        deriv = 4 * spins[i] * (norm(spins[i])^2 - S*(S+1))
        push!(derivs, deriv)
    end
    return derivs
end

function constrained_functional_derivs(spins, J, B::Vector{<:Real}, correlation_derivs, constraint, S, λ)
    derivs = mean_field_deriv(spins, J) .+ correlation_derivs(spins, J) .+ λ * constraint(spins, S)
    for i in eachindex(derivs)
        derivs[i] .+= B
    end
    return derivs
end

function constrained_functional_derivs(spins, J, B::Vector{<:AbstractVector}, correlation_derivs, constraint, S, λ)
    derivs = mean_field_deriv(spins, J) .+ correlation_derivs(spins, J) .+ λ * constraint(spins, S)
    derivs += B
    return derivs
end