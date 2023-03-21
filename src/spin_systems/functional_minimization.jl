include("functionals.jl")
include("functional_derivatives.jl")

function minimize(functional, derivative, spins_init, func_args, deriv_args, S=1/2, α=0.01, ϵ=1e-6, max_iter=100)
    E = functional(spins_init, func_args...)
    diff = 1.0
    iter = 1
    spins = spins_init
    while (diff > ϵ && iter <= max_iter)
        iter += 1
        derivs = derivative(spins, deriv_args...)
        spins_new = spins .- α * derivs
        spins_new = normalize_spin(spins_new, S)
        E_new = functional(spins_new, func_args...)
        diff = abs(E_new - E)
        spins = copy(spins_new)
        E = E_new
    end
    return E, spins, diff, iter
end


function example(N_spins::Int, S::Float64)
    spins_init = []
    for i in 1:N_spins
        push!(spins_init, rand(3))
    end
    spins_init = normalize_spin(spins_init, S)

    
    J = 1.0
    B = [rand(3) for i in 1:N_spins]
    α = 0.5
    λ = 1.0

    correlation = lsa_heisenberg_afm_1d
    corr_derivs = lsa_heisenberg_derivs
    constraint = spin_constraint_L2
    deriv_constraint = spin_constraint_L2_derivs

    func_args = (J, B, correlation, constraint, α, λ)
    deriv_Args = (J, B, corr_derivs, deriv_constraint, α, λ)

    functional = functional_with_constraint
    derivative = constrained_functional_derivs
        
    res = minimize(functional, derivative, spins_init, func_args, deriv_Args, S, α)
    return res
end