function two_qubit_entanglement(Z_expectations, k, l)
    N = length(Z_expectations)
    ent = 1 / (N-2)
    for j in [k,l]
        term = sum(Z_expectations) - Z_expectations[j]
        term += (N-3) * Z_expectations[j]
        ent *= sqrt(abs(term))
    end
    return ent
end