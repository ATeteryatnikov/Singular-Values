include("src/LogGenerator.jl")
include("lanczos.jl")
include("src/Matrices.jl")     
include("implicitQR.jl")     

const MIN_AVAILABLE_MARIX_SIZE = 6

global minMatrixSize
global prec
global limit

"""
Запускает алгоритм divide-and-conquer реализованный по статье:
"V. Rokhlin
A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices".

* Функция возвращает собственные числа поданной матрицы;
* Функция принимает:
	- alphas - главная диагональ;
	- bettas - диагональ соседняя с главной;
	- precB - критерий остановки метода бисекций, наибольшее отклонение характерстического полинома от нуля;
	- limitB - критерий остановки метода бисекций, наибольшее количество делений пополам;
	- minLength - минимальный размер матрицы при котором происходит разделение матрицы иначе
расчет собственных чисел производится неявным QR-алгоритмом.
"""
function divideAndConquer(alphas::AbstractArray, bettas::AbstractArray,
							precB, limitB, minLength)
	
	if (minLength < MIN_AVAILABLE_MARIX_SIZE)
		throw("Ошибка! Допустимый минимальный размер матрицы должен быть не меньше '6'!")
	end

	global minMatrixSize = minLength
	global prec = precB
	global limit = limitB
	
	return matrixPartsSVFinder(alphas, bettas, "")["T"]
end

"""
Вычисление собственных чисел плотной матрицы QR-алгоритмом со сдвигами Уилкинсона.

* Функция возвращает собственные числа.
* Функция принимает:
	- A - произвольная плотная матрица;
	- prec - максимальная норма соседней с главной диагонали матрицы.
"""
function denseQRAlgorithm(A, prec)
	if (size(A)[1] == 1)
		return A
	end
    A = copy(A)
    Q = Diagonal(ones(size(A)[1]))
    r = 1
    k = 1
    while (r > prec)
    	sub = A[end - 1:end,end - 1:end]
    	shift = WilkinsonShift(sub[1,1], sub[2,2], sub[1,2])
        (Q1, R) = qr(A - shift * I)
        A = R * Q1 + shift * I
        #Q = Q1'*Q
        r = norm(diag(A, -1))
    end
	return diag(A)
    #return [Q, A]
end

"""
Разделяет трехдиагональную матрицу в позиции [k,k] на две дочерние T*0 и T*1 

* Функция возвращает хэш-таблицу c элементом bettas[k] и 
с наборами сингулярных чисел подматриц с кодами:
T*binaryCode*0, T*binaryCode*0,
S*binaryCode*0, S*binaryCode*1,
R*binaryCode*0, R*binaryCode*1,
Q*binaryCode*0, Q*binaryCode*1.
* Функция принимает:
	- alphas - главная диагональ;
	- bettas - диагональ соседняя с главной;
	- binaryCode - бинарный код;
"""
function matrixHalfDivider(alphas, bettas, binaryCode)
	
	dict = Dict()
   	len = length(alphas)	

	k = div(len, 2)

	alphas[k] = alphas[k] - abs(bettas[k])
    alphas[k + 1] = alphas[k + 1] - abs(bettas[k])

    alphasT0 = alphas[1:k]
	bettasT0 = bettas[1:k-1]
    alphasT1 = alphas[k+1:end]
    bettasT1 = bettas[k+1:end]

    dictT0 = matrixPartsSVFinder(alphasT0, bettasT0, binaryCode * "0")
    dictT1 = matrixPartsSVFinder(alphasT1, bettasT1, binaryCode * "1")

	dict = merge(dictT1, dictT0)
    push!(dict, "b" => bettas[k])
	
	return dict
end

"""
Выполняет разделение матрицы через matrixHalfDivider,
получает собственные числа дочерних матриц и их частей,
по полученным наборам, через matrixPartsEvaluator, вычисляет 
собственные числа поданной на вход матрицы и её частей (T, S, R, Q).

* возвращает хэш-таблицу с наборами собственных чисел частей матрицы на уровне binaryCode с кодами:
 - T*binaryCode;
 - S*binaryCode;
 - R*binaryCode;
 - Q*binaryCode.
* Функция принимает:
	- alphas - главная диагональ;
	- bettas - диагональ соседняя с главной;
	- binaryCode - бинарный код;
"""
function matrixPartsSVFinder(alphas::AbstractArray, bettas::AbstractArray, binaryCode)
	dict = Dict()
	if (length(alphas)==1)
		push!(dict, "T" * binaryCode => spsvd(alphas, bettas))
		return dict
	end

	len = length(alphas)
	println("type: " * binaryCode)
	if (len < minMatrixSize)
		#spsvd = (a, b)->denseQRAlgorithm(toDense(a, b), prec)
		# Метод расчета собственных чисел для матриц меньше minMatrixSize
		# Могут возникать отрицательные собственные числа в ходе разделения матриц
		spsvd = (a, b)->RecursiveImplicitQR(copy(a), copy(b), prec, 200)
		
		push!(dict, "T" * binaryCode => spsvd(alphas, bettas))
		push!(dict, "S" * binaryCode => spsvd(alphas[1:end - 1], bettas[1:end - 1]))
		push!(dict, "R" * binaryCode => spsvd(alphas[2:end], bettas[2:end]))
		push!(dict, "Q" * binaryCode => spsvd(alphas[2:end - 1], bettas[2:end - 1]))
		
		return dict
	end

	submatrixDict = matrixHalfDivider(alphas, bettas, binaryCode)
	
	parentMatrixTSRQ = matrixPartsEvaluator(submatrixDict, binaryCode)
	
	return parentMatrixTSRQ
end

"""
Вычисляет собственные числа частей матрицы (T, S, R, Q) на уровне binaryCode.

* Функция принимает хэш-таблицу с наборами собственных чисел подматриц с кодами:
T*binaryCode*0, T*binaryCode*0,
S*binaryCode*0, S*binaryCode*1,
R*binaryCode*0, R*binaryCode*1,
Q*binaryCode*0, Q*binaryCode*1.
* Функция возвращает хэш-таблицу с наборами собственных чисел частей матрицы на уровне binaryCode с кодами:
 - T*binaryCode;
 - S*binaryCode;
 - R*binaryCode;
 - Q*binaryCode.
"""
function matrixPartsEvaluator(submatrixDict, binaryCode)

	bk = abs(submatrixDict["b"])

    # polyEval - функция вычисления значения полинома в точке, 
    # как произведения разностей собственных чисел и указанного x.
    # Функция принимает точку x и type - строка служащая ключем для хэш-таблицы,
    # если значение в хэш-таблицы не найдено - возвращает единицу т.к.
    # считаем, что значение характерстического полинома пустой матрицы единица
	function polyEval(x, type)
		if (haskey(submatrixDict, type))
			return prod(submatrixDict[type] .- x)
		else
			return 1
		end
	end
	evaluator = (x, X1, X2, X3, X4)->
		polyEval(x, X1) * polyEval(x, X2) + bk * (polyEval(x, X1) * polyEval(x, X4) + polyEval(x, X2) * polyEval(x, X3))

	Tevaluator = (x)->evaluator(x, "T" * binaryCode * "0", "T" * binaryCode * "1", "S" * binaryCode * "0", "R" * binaryCode * "1")
	Sevaluator = (x)->evaluator(x, "T" * binaryCode * "0", "S" * binaryCode * "1", "S" * binaryCode * "0", "Q" * binaryCode * "1")
	Revaluator = (x)->evaluator(x, "R" * binaryCode * "0", "T" * binaryCode * "1", "Q" * binaryCode * "0", "R" * binaryCode * "1")
	Qevaluator = (x)->evaluator(x, "R" * binaryCode * "0", "S" * binaryCode * "1", "Q" * binaryCode * "0", "Q" * binaryCode * "1")
	
	dict = Dict()
	
	hasZero = isa(findfirst("0", binaryCode), UnitRange)
	hasOne = isa(findfirst("1", binaryCode), UnitRange)

	# Условия необходимости вычисления части матрицы  
	# в источнике, Theorem 4.7.
	if (binaryCode == "")
		# evaluate only "T" matrix
	elseif (binaryCode == "0")
		singValsS = bisectBetween(submatrixDict["T" * binaryCode * "0"], submatrixDict["S" * binaryCode * "1"], bk, Sevaluator)
		push!(dict, "S" * binaryCode => singValsS)
	elseif (binaryCode == "1")
		singValsR = bisectBetween(submatrixDict["R" * binaryCode * "0"], submatrixDict["T" * binaryCode * "1"], bk, Revaluator)
		push!(dict, "R" * binaryCode => singValsR)
	else

		if (hasOne && hasZero)
			singValsQ = bisectBetween(submatrixDict["R" * binaryCode * "0"], submatrixDict["S" * binaryCode * "1"], bk, Qevaluator)
			push!(dict, "Q" * binaryCode => singValsQ)
		end
		if (hasOne)
			singValsR = bisectBetween(submatrixDict["R" * binaryCode * "0"], submatrixDict["T" * binaryCode * "1"], bk, Revaluator)
			push!(dict, "R" * binaryCode => singValsR)
		end
		if (hasZero)
			singValsS = bisectBetween(submatrixDict["T" * binaryCode * "0"], submatrixDict["S" * binaryCode * "1"], bk, Sevaluator)
			push!(dict, "S" * binaryCode => singValsS)
		end
	end
	singValsT = bisectBetween(submatrixDict["T" * binaryCode * "0"], submatrixDict["T" * binaryCode * "1"], bk, Tevaluator)
	push!(dict, "T" * binaryCode => singValsT)
	
	return dict
end

""" 
Поиск собственных чисел родительской матрицы в интервалах между собственными числами её дочерних матриц методом бисекций.

* Функция возвращает собственные числа родительской матрицы;
* Функция принимает:
	- eigvals1 - собственные числа первой матрицы;
	- eigvals2 - собственные числа второй матрицы;
	- bk - внедиагональный зануляемый при разделении родителькой матрицы элемент;
	- evaluator - функция вычисляющая значение характерстического полинома.
"""
function bisectBetween(eigvals1, eigvals2, bk, evaluator::Function)
	eigVals = []
	unionSV = sort(vcat(eigvals1, eigvals2))
	len = length(unionSV)
	for i in 1:1:len - 1
		push!(eigVals, bisection(evaluator, unionSV[i], unionSV[i + 1], prec, limit))
	end
	push!(eigVals, bisection(evaluator, unionSV[len], unionSV[len] + 2 * bk, prec, limit))

	return eigVals
end

"""
Вычисление корня характеристического полинома методом половинного 
деления на интервале.

* Функция возвращает корень характеристического полинома матрицы.
* Функция принимает:
	- rootFunc - функция вычисляющая значение характерстического полинома;
	- a, b - точки начала и конца интервала в котором происходит поиск корня;
	- prec - критерий остановки метода бисекций, наибольшее отклонение характерстического полинома от нуля;
	- limit - критерий остановки метода бисекций, наибольшее количество делений пополам.
"""
function bisection(rootFunc::Function, a, b, prec, limit)
	rc = 0
    c = 0
    ra = 0
    rb = 0
    a_changed = true
    b_changed = true
    k = 1
    while(k < limit)
        k = k + 1
        c = (a + b) / 2
        if (a_changed)
        	ra = rootFunc(a)
        	a_changed = false
        end
        if (b_changed)
        	rb = rootFunc(b)
        	b_changed = false
        end
        
        if (sign(ra) == sign(rb))
            return a
        end
        rc = rootFunc(c)

        if (sign(ra) == sign(rc))
        	a_changed = true
            a = c
        else
        	b_changed = true
            b = c
        end
        if (abs(rc) < prec)
            return c
        end
    end
    return c
end