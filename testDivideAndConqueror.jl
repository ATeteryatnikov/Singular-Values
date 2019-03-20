# Программа для тестирования алгоритма divide-and-conqueror
# из V. Rokhlin "A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices"
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("divideAndConqueror.jl") 
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 50
stepDim = 50
endDim = 100
# величина мантиссы BigFloat
startMantissa = 50
stepMantissa = 50
endMantissa = 100

startBisectionPrec = 10
stepBisectionPrec = 10
endBisectionPrec = 20

lenBlock = 10

# Директория с результатами теста
resultFolderName = "resultDivideAndConqueror"

# Создаем директорию в каталоге с программой с названием resultFolderName
try
	mkdir(string(pwd(), "\\", resultFolderName))
catch(err)
	println(err)
end

for n in startDim:stepDim:endDim
	for mantissa in startMantissa:stepMantissa:endMantissa
		for bisectionPrec in startBisectionPrec:stepBisectionPrec:endBisectionPrec

			println("Execute: dim=", n, " mantissa=", mantissa, "_bisectionPrec=1e-",bisectionPrec, "lenBlock=",lenBlock, "\n")

			setprecision(mantissa)
			godunovLists = lanczos(getGodunovMatrix(n)[1])
			alphas = godunovLists[1]
			bettas = godunovLists[2]
			fullMatr = toDense(alphas, bettas)
			
			DAC_svdvals, DACtimer = @timed svdValsFinder(copy(alphas), copy(bettas), 4, big(1)/big(10)^(bisectionPrec))
			
			println("divideAndConqueror complete!")

			Julia_svdvals, JuliaSVDTimer = @timed GenericSVD.generic_svdvals!(fullMatr)			
			
			difference = sort(Julia_svdvals,rev=true) - sort(abs.(DAC_svdvals),rev=true)
			normDifference = norm(difference)

			### Запись в файлы ###
			fileNameParameters = string(
				"_dim_", n, "_mantissa_", mantissa, "_bisectionPrec_1e-", bisectionPrec,"_lenBlock_", lenBlock, "_.txt")
			
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
	end
	end
end