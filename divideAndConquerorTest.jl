# Программа для тестирования алгоритма divide-and-conqueror
# из V. Rokhlin "A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices"
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("divideAndConqueror.jl") 
include("JacobiRotation.jl")
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 50
stepDim = 50
endDim = 50
# величина мантиссы BigFloat
startMantissa = 100
stepMantissa = 100
endMantissa = 100

startBisectionPrec = 1e-20
stepBisectionPrec = 1
endBisectionPrec = 1e-20

lenBlock = 4

# Директория с результатами теста
resultFolderName = ""

# Если не указан resultFolderName, то директория будет называться "resultJacobiRotationTests"
if (resultFolderName=="")
	
	resultFolderName = "resultDivideAndConquerorTest"
end
# Создаем директорию в каталоге с программой с названием resultFolderName
try
	mkdir(string(pwd(), "\\", resultFolderName))
catch(err)
	println(err)
end

for n in startDim:stepDim:endDim
	for mantissa in startMantissa:stepMantissa:endMantissa
		for bisectionPrec in startBisectionPrec:stepBisectionPrec:endBisectionPrec

			println("Execute: dim=", n, " mantissa=", mantissa, "_bisectionPrec=",bisectionPrec, "\n")

			setprecision(mantissa)
			godunovLists = lanczos(getGodunovMatrix(n)[1])
			alphas = godunovLists[1]
			bettas = godunovLists[2]
			fullMatr = toDense(alphas, bettas)
			
			DAC_svdvals, DACtimer = @timed svdValsFinder(copy(alphas), copy(bettas), 4, bisectionPrec)
			
			Julia_svdvals, JuliaSVDTimer = @timed GenericSVD.generic_svdvals!(fullMatr)			
			
			difference = sort(Julia_svdvals,rev=true) - sort(abs.(DAC_svdvals),rev=true)
			normDifference = norm(difference)

			### Запись в файлы ###
			# Разность сингулярных чисел полученных методом divide-and-conqueror и командой svd.
			pathDiffSingVals = string(resultFolderName,"/DiffSingVals__dim_", n, "_mantissa_", mantissa, "_bisectionPrec_", bisectionPrec,"_.txt")
			fileDifferenceSingVals = open(pathDiffSingVals, "w")
			writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],","))
			close(fileDifferenceSingVals)

			# Норма разности сингулярных чисел полученных методом divide-and-conqueror и командой svd.
			pathNormDiffSingVals = string(resultFolderName,"/NormDiffSingVals__dim_", n, "_mantissa_", mantissa, "_bisectionPrec_", bisectionPrec, "_.txt")
			fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
			write(fileNormDifferenceSingVals, string(normDifference))
			close(fileNormDifferenceSingVals)

			# Время работы метода divide-and-conqueror
			pathDACTime = string(resultFolderName,"/DACTime__dim_", n, "_mantissa_", mantissa, "_bisectionPrec_", bisectionPrec, "_.txt")
			fileDACTime = open(pathDACTime, "w")
			writedlm(fileDACTime, DACtimer)
			close(fileDACTime)

			# Время работы команды svd
			pathSVDTime = string(resultFolderName,"/JuliaSVDTime__dim_", n, "_mantissa_", mantissa, "_bisectionPrec_", bisectionPrec, "_.txt")
			fileSVDTime = open(pathSVDTime, "w")
			writedlm(fileSVDTime, JuliaSVDTimer)
			close(fileSVDTime)

			pathInfoLog = string(resultFolderName,"/InfoLog__dim_", n, "_mantissa_", mantissa, "_bisectionPrec_", bisectionPrec, "_.txt")
			fileInfoLog = open(pathInfoLog, "w")
			write(fileInfoLog, string("Мантисса : ", mantissa, "\n"))
			write(fileInfoLog, string("Размерность СЛАУ Годунова : ", n, "\n"))
			write(fileInfoLog, string("Время работы команды svd : ", JuliaSVDTimer, "\n"))
			write(fileInfoLog, string("Время работы метода divide-and-conqueror : ", DACtimer, "\n"))
			write(fileInfoLog, string("Точность вычисления в методе бисекций : ", bisectionPrec, "\n"))
			write(fileInfoLog, string("Норма разности сингулярных чисел полученных методом divide-and-conqueror и командой svd : ", normDifference, "\n"))
			close(fileInfoLog)
	end
	end
end