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



