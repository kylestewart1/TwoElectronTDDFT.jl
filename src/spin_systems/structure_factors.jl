function correlation(Ψ, paulis, indices=(1,2), direction="z")
    i, j = indices
    op = paulis[direction][i]*paulis[direction][j]
    return Ψ' * op * Ψ
end


function  correlation_function(Ψ, paulis, n_qubits, direction="z")
    corr = []
    for i in 2:n_qubits
        push!(corr, correlation(Ψ,paulis,(i,1),direction))
    end
    return corr
end

function structure_factor(p, corr)
    struc = zeros(ComplexF64, size(p))
    for l in eachindex(corr)
        struc .+= corr .* exp.(im*p*l)
    end
    return struc
end
