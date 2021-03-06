# Программа для тестирования неявного QR-алгоритма.
# Watkins D. S., Elsner L. Chasing algorithms for the eigenvalue problem 
# //SIAM journal on matrix analysis and 
# applications. – 1991. – Т. 12. – №. 2. – С. 374-384
# Матрица для тестов - левая часть СЛАУ Годунова преобразованная к трехдиагональному виду методом Ланцоша.

include("lanczos.jl")
include("implicitQR.jl")
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles
# размерность СЛАУ Годунова
startDim = 50
stepDim = 200
endDim = 50
# величина мантиссы BigFloat
startMantissa = 100
stepMantissa = 50
endMantissa = 100
# Пороговое значение зануления внедиагональных элементов
startPrecision = 40
stepPrecision = 40
endPrecision = 80
# Количество итераций неявного QR-алгоритма между проведением зануления внедиагональных элементов
startNumIteration = 200
stepNumIteration = 200
endNumIteration = 200

# Директория с результатами теста
resultFolderName = ""

# Если не указан resultFolderName, то директория будет называться "resultJacobiRotationTests"
resultFolderName = "resultImplicitQR11"

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
				implQRVals, implQRTimer = @timed RecursiveImplicitQR(copy(alphas), copy(bettas), big(1)/big(10)^(prec), numIteration)
	
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