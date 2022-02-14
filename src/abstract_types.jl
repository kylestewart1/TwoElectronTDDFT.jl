import Base: +, *, isapprox

# numerical type for elements of all quantum objects
const Scalar = Union{Float64, ComplexF64} 

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


