# Программа для тестирования классического метода Якоби.
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("Lanczos.jl")
include("ImplicitQR.jl")
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 100
stepDim = 50
endDim = 300
# величина мантиссы BigFloat
startMantissa = 100
stepMantissa = 100
endMantissa = 100
# Пороговое значение зануления внедиагональных элементов
startPrecision = 80
stepPrecision = 80
endPrecision = 80
# Количество итераций неявного QR-алгоритма между проведением зануления внедиагональных элементов
startNumIteration = 50
stepNumIteration = 50
endNumIteration = 50

# Директория с результатами теста
resultFolderName = ""

# Если не указан resultFolderName, то директория будет называться "resultJacobiRotationTests"
if (resultFolderName=="")
	# resultFolderName = "resultJacobiRotationTests"
	resultFolderName = "resultImplicitQR"
end
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
		println("Расчет svdvals")
		originVals, denseTimer = @timed svdvals(fullMatr)
		println("Расчет svdvals окончен ",denseTimer," сек.")

		for prec in startPrecision:stepPrecision:endPrecision
			for numIteration in startNumIteration:stepNumIteration:endNumIteration

				println("\nExecute: dim=", n, " mantissa=", mantissa, " prec=1e-", prec, " numIteration=", numIteration)

				
				print("Расчет implQR")
				implQRVals, implQRTimer = @timed RecursiveImplicitQR(copy(alphas), copy(bettas), big(1)/big(10)^(prec), numIteration, 10000)
				implQRVals = BigFloat.(implQRVals[1])
				print("Расчет implQR окончен ",implQRTimer," сек.")
				
				
				difference = sort(originVals,rev=true) - sort(abs.(implQRVals),rev=true)
				normDifference = norm(difference)
			
				### Запись в файлы ###
				# Разность сингулярных чисел полученных классическим методом Якоби и командой svd.
				pathDiffSingVals = string(resultFolderName,"/DiffSingVals__dim_", n, "_mantissa_", mantissa, "_prec_1e", -prec,"_numIteration_", numIteration,"_.txt")
				fileDifferenceSingVals = open(pathDiffSingVals, "w")
				writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],","))
				close(fileDifferenceSingVals)

				# Норма разности сингулярных чисел полученных классическим методом Якоби и командой svd.
				pathNormDiffSingVals = string(resultFolderName,"/NormDiffSingVals__dim_", n, "_mantissa_", mantissa, "_prec_1e", -prec,"_numIteration_", numIteration, "_.txt")
				fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
				write(fileNormDifferenceSingVals, string(normDifference))
				close(fileNormDifferenceSingVals)

				# Время работы классического метода Якоби
				pathImplQRTime   = string(resultFolderName,"/ImplQRTime__dim_", n, "_mantissa_", mantissa, "_prec_1e", -prec,"_numIteration_", numIteration, "_.txt")
				fileImplQRTime = open(pathImplQRTime, "w")
				write(fileImplQRTime, string(implQRTimer))
				close(fileImplQRTime)

				# Время работы команды svd
				pathDenseTime = string(resultFolderName,"/DenseTime__dim_", n, "_mantissa_", mantissa, "_prec_1e", -prec,"_numIteration_", numIteration, "_.txt")
				fileDenseTime = open(pathDenseTime, "w")
				writedlm(fileDenseTime, denseTimer)
				close(fileDenseTime)


				pathInfoLog = string(resultFolderName,"/InfoLog__dim_", n, "_mantissa_", mantissa, "_prec_1e", -prec,"_numIteration_", numIteration, "_.txt")
				fileInfoLog = open(pathInfoLog, "w")
				write(fileInfoLog, string("Мантисса : ", mantissa, "\n"))
				write(fileInfoLog, string("Размерность СЛАУ Годунова : ", n, "\n"))
				write(fileInfoLog, string("Точность : ", prec, "\n"))
				write(fileInfoLog, string("Итераций между дефляцией : ", numIteration, "\n"))
				write(fileInfoLog, string("Время работы команды svd : ", denseTimer, "\n"))
				write(fileInfoLog, string("Время работы неявного QR-алгоритма : ", implQRTimer, "\n"))
				write(fileInfoLog, string("Норма разности сингулярных чисел полученных неявным QR-алгоритмом и командой svd : ", normDifference, "\n"))
				close(fileInfoLog)
			end
		end
	end
end