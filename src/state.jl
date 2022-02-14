

struct Ket <: AbstractState
    data::Vector{Scalar}
end

struct Bra <: AbstractState
    data::Vector{Scalar}
end

adjoint(ψ::Ket) = Bra(conj(data(ψ)))
adjoint(ϕ::Bra) = Ket(conj(data(ϕ)))

*(ϕ::Bra, ψ::Ket) = check_dimensions_match(ϕ, ψ) && data(ϕ) ⋅ data(ψ)
*(ψ::Ket, ϕ::Bra) = check_dimensions_match(ψ, ϕ) && Operator(data(ψ) * data(ϕ)')