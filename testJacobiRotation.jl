# Программа для тестирования классического метода Якоби.
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("Lanczos.jl")
include("JacobiRotation.jl")
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 300
stepDim = 300
endDim = 300
# величина мантиссы BigFloat
startMantissa = 100
stepMantissa = 100
endMantissa = 100
# Кол-во итераций в алгоритме
startNumIterations = 20000
stepNumIterations = 1000
endNumIterations = 20000
# Директория с результатами теста
resultFolderName = ""
# в алгоритме в диагонали bettas ищется первый элемент болше чем maxBettasLimit с которого начинается итерация алгоритма
# если все элементы из bettas не превышают значения maxBettasLimit, то алгоритм прекращает работу
maxBettasLimit = 1e-40

# Если не указан resultFolderName, то директория будет называться "resultJacobiRotationTests"
if (resultFolderName=="")
	# resultFolderName = "resultJacobiRotationTests"
	resultFolderName = "resultJacobiRotationWithoutShiftTests"
end
# Создаем директорию в каталоге с программой с названием resultFolderName
try
	mkdir(string(pwd(), "\\", resultFolderName))
catch(err)
	println(err)
end

for n in startDim:stepDim:endDim
	for mantissa in startMantissa:stepMantissa:endMantissa
		for numRepeated in startNumIterations:stepNumIterations:endNumIterations

			println("Execute: dim=", n, " mantissa=", mantissa, " rep=", numRepeated)

			setprecision(mantissa)
			godunovLists = lanczos(getGodunovMatrix(n)[1])
			alphas = godunovLists[1]
			bettas = godunovLists[2]
			fullMatr = toDense(alphas, bettas)
			
			# JacobiRotation возвращает два массива (диагонали) и количество выполненых итераций [alphas, bettas, k]
			JRRes, timerJR = @timed JacobiRotation(copy(alphas), copy(bettas), numRepeated)
			#JRRes, timerJR = @timed JacobiRotationModification(copy(alphas), copy(bettas), numRepeated, maxBettasLimit)
			#JRRes, timerJR = @timed JacobiRotationWithShift(copy(alphas), copy(bettas), numRepeated, maxBettasLimit)
			
			originVals, denseTimer = @timed svdvals(fullMatr)
			
			difference = sort(originVals,rev=true) - sort(abs.(JRRes[1]),rev=true)
			normDifference = norm(difference)
			maxBetta = maximum(JRRes[2])

			maxBettasLimitComplete = false
			#if (numRepeated == JRRes[3])
			#	maxBettasLimitComplete = true
			#end

			### Запись в файлы ###
			# Разность сингулярных чисел полученных классическим методом Якоби и командой svd.
			pathDiffSingVals = string(resultFolderName,"/DiffSingVals__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileDifferenceSingVals = open(pathDiffSingVals, "w")
			writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],","))
			close(fileDifferenceSingVals)

			# Норма разности сингулярных чисел полученных классическим методом Якоби и командой svd.
			pathNormDiffSingVals = string(resultFolderName,"/NormDiffSingVals__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
			write(fileNormDifferenceSingVals, string(normDifference))
			close(fileNormDifferenceSingVals)

			# Время работы классического метода Якоби
			pathJRTime = string(resultFolderName,"/JRTime__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileJRTime = open(pathJRTime, "w")
			writedlm(fileJRTime, timerJR)
			close(fileJRTime)

			# Время работы команды svd
			pathDenseTime = string(resultFolderName,"/DenseTime__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileDenseTime = open(pathDenseTime, "w")
			writedlm(fileDenseTime, denseTimer)
			close(fileDenseTime)

			# Вектор bettas
			pathBettas = string(resultFolderName,"/Bettas__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileBettas = open(pathBettas, "w")
			writedlm(fileBettas, split(string(JRRes[2])[10:end-1],","))
			close(fileBettas)

			# максимальный эллемент из bettas
			pathMaxBetta = string(resultFolderName,"/MaxBetta__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileMaxBetta = open(pathMaxBetta, "w")
			write(fileMaxBetta, string(maxBetta))
			close(fileMaxBetta)

			pathInfoLog = string(resultFolderName,"/InfoLog__dim_", n, "_mantissa_", mantissa, "_rep_", numRepeated, "_.txt")
			fileInfoLog = open(pathInfoLog, "w")
			write(fileInfoLog, string("Мантисса : ", mantissa, "\n"))
			write(fileInfoLog, string("Размерность СЛАУ Годунова : ", n, "\n"))
			write(fileInfoLog, string("Заданное количество итераций : ", numRepeated, "\n"))
			#write(fileInfoLog, string("Итераций выполнено : ", JRRes[3], "\n"))
			write(fileInfoLog, string("Время работы команды svd : ", denseTimer, "\n"))
			write(fileInfoLog, string("Время работы классического метода Якоби : ", timerJR, "\n"))
			write(fileInfoLog, string("максимальный эллемент из соседней с главной диагонали : ", maxBetta, "\n"))
			write(fileInfoLog, string("Норма разности сингулярных чисел полученных классическим методом Якоби и командой svd : ", normDifference, "\n"))
			write(fileInfoLog, string("Было выполнено все заданное количество итераций : ", maxBettasLimitComplete, "\n"))
			close(fileInfoLog)
		end
	end
end