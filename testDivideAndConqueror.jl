# Программа для тестирования алгоритма divide-and-conqueror
# из V. Rokhlin "A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices"
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("divideAndConqueror.jl") 
include("src/Matrices.jl")
using GenericSVD
using GenericLinearAlgebra
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 300
stepDim = 300
endDim = 300
# величина мантиссы BigFloat
startMantissa = 300
stepMantissa = 200
endMantissa = 300

startBisectionPrec = 40
stepBisectionPrec = 40
endBisectionPrec = 40

lenBlock = 50

# Директория с результатами теста
resultFolderName = "resultDivideAndConquerorTEST"

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
		
		for bisectionPrec in startBisectionPrec:stepBisectionPrec:endBisectionPrec

			println("Execute: dim=", n, " mantissa=", mantissa, "_bisectionPrec=1e-",bisectionPrec, "_lenBlock=",lenBlock, "\n")

			
			
			DAC_svdvals, DACtimer = @timed svdValsFinder(copy(alphas), copy(bettas), lenBlock, big(1)/big(10)^(bisectionPrec))
			
			println("divideAndConqueror complete!")

			Julia_svdvals, JuliaSVDTimer = @timed GenericSVD.generic_svdvals!(fullMatr)			
			
			difference = sort(Julia_svdvals,rev=true) - sort(abs.(DAC_svdvals),rev=true)
			normDifference = norm(difference)

			### Запись в файлы ###
			fileNameParameters = string(
				"_dim_", n, "_mantissa_", mantissa, "_bisectionPrec_1e-", bisectionPrec,"_lenBlock_", lenBlock, "_.txt")
			
			# Сингулярные числа найденные методом divide-and-conquer
			pathDiffSingVals = string(resultFolderName,"/SingVals_", fileNameParameters)
			fileDifferenceSingVals = open(pathDiffSingVals, "w")
			writedlm(fileDifferenceSingVals, split(string(DAC_svdvals)[10:end-1],","))
			close(fileDifferenceSingVals)

			# Разность сингулярных чисел полученных методом divide-and-conqueror и командой svd.
			pathDiffSingVals = string(resultFolderName,"/DiffSingVals_", fileNameParameters)
			fileDifferenceSingVals = open(pathDiffSingVals, "w")
			writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],","))
			close(fileDifferenceSingVals)

			# Норма разности сингулярных чисел полученных методом divide-and-conqueror и командой svd.
			pathNormDiffSingVals = string(resultFolderName,"/NormDiffSingVals_", fileNameParameters)
			fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
			write(fileNormDifferenceSingVals, string(normDifference))
			close(fileNormDifferenceSingVals)

			# Время работы метода divide-and-conqueror
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
			write(fileInfoLog, string("Время работы метода divide-and-conqueror : ", DACtimer, "\n"))
			write(fileInfoLog, string("Точность вычисления в методе бисекций : 1e-", bisectionPrec, "\n"))
			write(fileInfoLog, string("Норма разности сингулярных чисел полученных методом divide-and-conqueror и командой svd : ", normDifference, "\n"))
			close(fileInfoLog)

			println("Норма разности сингулярных чисел полученных методом svdvals и методом divide-and-conquer : ",normDifference)
	end
	end
end