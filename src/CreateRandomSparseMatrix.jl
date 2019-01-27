function createRandomSparseMatrix(m, n, countNotNullValues, useBigFloat = false)
    I1 = zeros(countNotNullValues)
    I1[end] = m
    I1[1:end-1] = rand(1:m, countNotNullValues-1)

    I2 = zeros(countNotNullValues)
    I2[end] = n
    I2[1:end-1] = rand(1:n, countNotNullValues-1)

    if useBigFloat
        return sparse(I1, I2, rand(big(-50.0):big(50.0), countNotNullValues))
    end
    return sparse(I1, I2, rand(1:10000, countNotNullValues))

end
