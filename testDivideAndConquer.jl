# Программа для тестирования алгоритма divide-and-conquer
# из V. Rokhlin "A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices"
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("divideAndConquer.jl") 
include("src/Matrices.jl")
using GenericSVD
using GenericLinearAlgebra
using DelimitedFiles

# размерность СЛАУ Годунова
startDim = 300
stepDim = 200
endDim = 300
# величина мантиссы BigFloat
startMantissa = 100
stepMantissa = 100
endMantissa = 100
# Наибольшее отклонение характеристического полинома от нуля в методе бисекций.
# Если значение полинома меньше заданного числа,
# считаем, что корень найден;
startBisectionPrec = 40
stepBisectionPrec = 40
endBisectionPrec = 80
# максимальная размерность матриц в рекурсивном разделении матриц, 
# при которой расчет значения характеристического полинома происходит через
# функцию Julia "det".
startLenBlock = 6
stepLenBlock = 6
endLenBlock = 6


# Директория с результатами теста
resultFolderName = "resultDivideAndConquerExperimental"

# Создаем директорию в каталоге с программой с названием resultFolderName
try
	mkdir(string(pwd(), "\\", resultFolderName))
catch(err)
	println(err)
end

for n in startDim:stepDim:endDim
	for mantissa in startMantissa:stepMantissa:endMantissa
		
		setprecision(mantissa)
		godunovLists = lanczos(getGodunovMatrix(n)[1])
		alphas = godunovLists[1]
		bettas = godunovLists[2]
		fullMatr = toDense(alphas, bettas)

		Julia_svdvals, JuliaSVDTimer = @timed GenericSVD.generic_svdvals!(fullMatr)

		for bisectionPrec in startBisectionPrec:stepBisectionPrec:endBisectionPrec
			for lenBlock in startLenBlock:stepLenBlock:endLenBlock
			
				println("Execute: dim=", n, " mantissa=", mantissa, "_bisectionPrec=1e-",bisectionPrec, "_lenBlock=", lenBlock, "\n")

				DAC_svdvals, DACtimer = @timed divideAndConquer(copy(alphas), copy(bettas), big(1)/big(10)^(bisectionPrec), 500, lenBlock)
			
				println("divideAndConquer complete!")

						
				
				difference = sort(Julia_svdvals,rev=true) - sort(DAC_svdvals,rev=true)
				normDifference = norm(difference)

				### Запись в файлы ###
				fileNameParameters = string(
					"_dim_", n, "_mantissa_", mantissa, "_bisectionPrec_1e-", bisectionPrec,"_lenBlock_", lenBlock, "_.txt")
				
				# Сингулярные числа найденные методом divide-and-conquer
				pathDiffSingVals = string(resultFolderName,"/SingVals_", fileNameParameters)
				fileDifferenceSingVals = open(pathDiffSingVals, "w")
				writedlm(fileDifferenceSingVals, split(string(big.(DAC_svdvals))[10:end-1],", "))
				close(fileDifferenceSingVals)

				# Разность сингулярных чисел полученных методом divide-and-conquer и командой svd.
				pathDiffSingVals = string(resultFolderName,"/DiffSingVals_", fileNameParameters)
				fileDifferenceSingVals = open(pathDiffSingVals, "w")
				writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],", "))
				close(fileDifferenceSingVals)

				# Норма разности сингулярных чисел полученных методом divide-and-conquer и командой svd.
				pathNormDiffSingVals = string(resultFolderName,"/NormDiffSingVals_", fileNameParameters)
				fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
				write(fileNormDifferenceSingVals, string(normDifference))
				close(fileNormDifferenceSingVals)

				# Время работы метода divide-and-conquer
				pathDACTime = string(resultFolderName,"/DACTime_", fileNameParameters)
				fileDACTime = open(pathDACTime, "w")
				writedlm(fileDACTime, DACtimer)
				close(fileDACTime)

				# Время работы команды svd
				pathSVDTime = string(resultFolderName,"/JuliaSVDTime_", fileNameParameters)
				fileSVDTime = open(pathSVDTime, "w")
				writedlm(fileSVDTime, JuliaSVDTimer)
				close(fileSVDTime)

				pathInfoLog = string(resultFolderName,"/InfoLog_", fileNameParameters)
				fileInfoLog = open(pathInfoLog, "w")
				write(fileInfoLog, string("Мантисса : ", mantissa, "\n"))
				write(fileInfoLog, string("Размерность СЛАУ Годунова : ", n, "\n"))
				write(fileInfoLog, string("Время работы команды svd : ", JuliaSVDTimer, "\n"))
				write(fileInfoLog, string("Время работы метода divide-and-conquer : ", DACtimer, "\n"))
				write(fileInfoLog, string("Точность вычисления в методе бисекций : 1e-", bisectionPrec, "\n"))
				write(fileInfoLog, string("Норма разности сингулярных чисел полученных методом divide-and-conquer и командой svd : ", normDifference, "\n"))
				close(fileInfoLog)

				println("Норма разности сингулярных чисел полученных методом svdvals и методом divide-and-conquer : ",normDifference)
			end
		end
	end
end