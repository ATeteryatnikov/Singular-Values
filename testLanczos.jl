# Программа для тестирования метода Ланцоша

try
    include("lanczos.jl")
    include("src/Matrices.jl")
    using GenericSVD
catch ex
    print("Что-то пошло не так: ", ex)
end

# Мантисса чисел BigFloat
bigFloatMantissa = 100
# Размерность СЛАУ Годунова
n = 80

setprecision(bigFloatMantissa)

A = getGodunovMatrix(n)[1]

listsAlphAndBetta = lanczos(A)
denseTridiagonalMatrix = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])

###     Нахождение сигулярных чисел исходной и треугольной матриц
singValsOriginalMatrix = svdvals(A'*A)
singValsTridiagonalMatrix = svdvals(denseTridiagonalMatrix)

diffSingVals = sqrt.(singValsOriginalMatrix) - sqrt.(singValsTridiagonalMatrix)
normDiff = norm(diffSingVals)
print(Float64.(normDiff))

###     Запись в файлы
path = string("resultLanczosTests/diffSingVals_dim_", n, "_matissa_", bigFloatMantissa, ".txt")
fileDifferenceSingVals = open(path, "w")
# Записываем результат в формате Float64
writedlm(fileDifferenceSingVals, Float64.(diffSingVals))
close(fileDifferenceSingVals)

# сингулярные числа трехдиагональной матрицы
path = string("resultLanczosTests/singValsTridiagonalMatrix_dim_", n, "_matissa_", bigFloatMantissa, ".txt")
fileSingVals = open(path, "w")
# Записываем результат в формате Float64
writedlm(fileSingVals, Float64.(sqrt.(singValsTridiagonalMatrix)))
close(fileSingVals)

# сингулярные числа исходной матрицы
path = string("resultLanczosTests/singValsOriginMatrix_dim_", n, "_matissa_", bigFloatMantissa, ".txt")
fileSingValsOrigin = open(path, "w")
# Записываем результат в формате Float64
writedlm(fileSingValsOrigin, Float64.(sqrt.(singValsOriginalMatrix)))
close(fileSingValsOrigin)

# норма разности сингулярных чисел трехдиагональной матрицы и исходной матрицы домноженной слева на транспонированную
path = string("resultLanczosTests/normDiffSingVals_dim_", n, "_matissa_", bigFloatMantissa, ".txt")
fileNormDiffSingVals = open(path, "w")
# Записываем результат в формате Float64
write(fileNormDiffSingVals, string(normDiff))
close(fileNormDiffSingVals)
