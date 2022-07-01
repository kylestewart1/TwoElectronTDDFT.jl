import Base: +, *, isapprox


abstract type AbstractQuantumObject end
abstract type AbstractQuantumMatrix <: AbstractQuantumObject end
abstract type AbstractQuantumVector <: AbstractQuantumMatrix end
abstract type AbstractPropagator end


data(A::AbstractQuantumObject) = A.data
size(A::AbstractQuantumObject) = size(A.data)
getindex(A::AbstractQuantumObject, i::Integer) = A.data[i]
adjoint(A::AbstractQuantumObject) = error("Adjoint is not defined for abstract type", AbstractQuantumObject)

isapprox(A::AbstractQuantumMatrix, obj2::AbstractQuantumMatrix) = size(obj1) == size(obj2) && obj1.data ≈ obj2.data
adjoint(op::AbstractQuantumMatrix) = adjoint(op.data)
ishermitian(op::AbstractQuantumMatrix) = ishermitian(op.data)
isunitary(op::AbstractQuantumMatrix) = isunitary(op.data)

norm(ψ::AbstractQuantumMatrix) = norm(ψ.data)

+(obj1::AbstractQuantumObject, obj2::AbstractQuantumObject) = error("Addition not defined for abstract type", typeof(obj1))
*(obj1::AbstractQuantumObject, obj2::AbstractQuantumObject) = error("Multiplication not defined for abstract type", typeof(obj1))

function check_dimensions_match(obj1::AbstractQuantumObject, obj2::AbstractQuantumObject)
    size(obj1)==size(obj2) || throw(DimensionMismatch("Dimensions of ", obj1, " and ", obj2, " do not match."))
end

function mul!(C::AbstractQuantumMatrix, A::AbstractQuantumMatrix, B::AbstractQuantumMatrix, α::Number, β::Number)
    mul!(data(C), data(A), data(B), α, β)
end


function update!(A::AbstractQuantumObject, data)
    for i in eachindex(data)
        A.data[i] = data[i]
    end
end


"""
Operator{T <: AbstractMatrix} <: AbstractQuantumMatrix

Quantum operator type. Elements are stored in field `data` as a matrix of type T.
"""
struct Operator{T <: AbstractMatrix} <: AbstractQuantumMatrix
    data::T
    function Operator{T}(A::AbstractArray) where {T}
        new{T}(T(A))
    end    
end

Operator(A::T) where {T <: AbstractMatrix} = Operator{T}(A)
Operator(A::T) where {T <: AbstractVector} = Operator(Diagonal(A))
# `zeros`-like constructor for NxN operator
Operator(N::Integer) = Operator(zeros(N, N))

# type aliases for convenience and clarity
const Hamiltonian{T} = Operator{Hermitian{T}}
const Potential = Operator{Diagonal{Float64}}
const KineticEnergy{T} = Operator{Tridiagonal{T}}

"""
Ket{T <: AbstractVector} <: AbstractState

An element of the Hilbert space. To be thought of as a column vector,
but elements `data` can be any vector type.
"""
struct Ket{T <: AbstractVector} <: AbstractQuantumVector
    data::T
end

"""
Bra{T <: AbstractVector} <: AbstractState

An element of the dual space. To be thought of as a row vector,
but elements `data` can be any vector type.
"""
struct Bra{T <: AbstractVector} <: AbstractQuantumVector
    data::T
end


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


adjoint(ψ::Ket) = Bra(conj(data(ψ)))
adjoint(ϕ::Bra) = Ket(conj(data(ϕ)))

*(ϕ::Bra, ψ::Ket) = check_dimensions_match(ϕ, ψ) && data(ϕ) ⋅ data(ψ)
*(ψ::Ket, ϕ::Bra) = check_dimensions_match(ψ, ϕ) && Operator(data(ψ) * data(ϕ)')


