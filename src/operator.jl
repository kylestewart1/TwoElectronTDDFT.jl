struct Operator{T <: AbstractMatrix{Scalar}} <: AbstractQuantumMatrix
    data::T
    function Operator{T}(A::AbstractMatrix)
        new{T}(A)
    end    
end

const Hamiltonian = Operator{Hermitian}
const Potential = Operator{Diagonal}
const KineticEnergy = Operator{Tridiagonal}

Operator(A::T) where {T <: AbstractMatrix{Scalar}} = Operator{T}(A)
Operator(N::Integer) = Operator(zeros(Scalar, N, N))


function (A::Operator)(ψ::Ket)
    check_dimensions_match(A,ψ)
    Ket(data(A) * data(ψ))
end

function (A::Operator)(ϕ::Bra) 
    check_dimensions_match(A, ϕ)
    Bra(adjoint(data(A)) * data(ϕ))
end

function *(A::Operator, B::Operator)
    check_dimensions_match(A, B)
    Operator(data(A) * data(B))
end

*(A::Operator, ψ::Ket) = (A)(ψ)
*(ϕ::Bra, A::Operator) = (A)(ϕ)

function mul!(C::AbstractQuantumMatrix, A::AbstractQuantumMatrix, B::AbstractQuantumMatrix, α::Scalar, β::Scalar)
    mul!(data(C), data(A), data(B), α, β)
end




