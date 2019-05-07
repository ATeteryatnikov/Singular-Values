include("src/LogGenerator.jl")
include("lanczos.jl")
include("src/Matrices.jl")     
include("ImplicitQR.jl")     

"""
Вычисление собственных чисел плотной матрицы QR-алгоритмом со сдвигами Уилкинсона
"""
function denseQRAlgorithm(A, prec)
	if (size(A)[1]==1)
		return A
	end
    A = copy(A)
    Q = Diagonal(ones(size(A)[1]))
    r=1
    while (r>prec)
    	sub = A[end-1:end,end-1:end]
    	shift = WilkinsonShift(sub[1,1],sub[2,2],sub[1,2])
        (Q1, R) = qr(A - shift*I)
        A = R * Q1 + shift*I
        #Q = Q1'*Q
        r = norm(diag(A,-1))
    end
return diag(A)
    #return [Q, A]
end

minLen = 6

# Разделяет матрицу на две дочерние T0 и T1
# Отдает сингулярные числа всех подматриц 
function T0T1divider(alphas, bettas, binaryCode, precB, limitB)
	
	println("T0T1divider: "*binaryCode)

	alpha = copy(alphas)
	betta = copy(bettas)

	dict = Dict()
	len = length(alpha)	

	k = div(len, 2)

	alpha[k] = alpha[k] - abs(betta[k])
    alpha[k + 1] = alpha[k + 1] - abs(betta[k])

    alphasT0 = copy(alpha[1:k])
    bettasT0 = copy(betta[1:k-1])
    alphasT1 = copy(alpha[k+1:end])
    bettasT1 = copy(betta[k+1:end])

    dictT0 = SingValsFinder(alphasT0, bettasT0, binaryCode*"0", precB, limitB)
    dictT1 = SingValsFinder(alphasT1, bettasT1, binaryCode*"1", precB, limitB)

    dict = merge(dictT1, dictT0)
	push!(dict, "b" => betta[k])
    return dict
end	

function SingValsFinder(alphas, bettas, binaryCode, precB, limitB)
	len = length(alphas)
	println("type: "*binaryCode)
	if (len < minLen)
		spsvd = (a, b) -> denseQRAlgorithm(toDense(a, b), precB)
		dict = Dict()
		push!(dict, "T"*binaryCode => spsvd(alphas, bettas))
		push!(dict, "S"*binaryCode => spsvd(alphas[1:end-1], bettas[1:end-1]))
		push!(dict, "R"*binaryCode => spsvd(alphas[2:end], bettas[2:end]))
		push!(dict, "Q"*binaryCode => spsvd(alphas[2:end-1], bettas[2:end-1]))
		return dict
	end

	submatrixDict = T0T1divider(alphas, bettas, binaryCode, precB, limitB)
	
	parentMatrixTSRQ = parentMatrixEvaluator(submatrixDict, binaryCode, precB, limitB)
	
	return parentMatrixTSRQ
end

function parentMatrixEvaluator(submatrixDict, binaryCode, precB, limitB)

	bk = abs(submatrixDict["b"])

	# функция вычисления полинома в точке по собственным числам дочерних матриц
	function polyEval(x, type)
		if (haskey(submatrixDict, type))
			return prod(submatrixDict[type] .- x)
		else
			return 1
		end
	end
	evaluator = (x, X1, X2, X3, X4) -> polyEval(x, X1)*polyEval(x, X2) + bk*(polyEval(x, X1)*polyEval(x, X4) + polyEval(x, X2)*polyEval(x, X3))

	Tevaluator = (x) -> evaluator(x, "T"*binaryCode*"0", "T"*binaryCode*"1", "S"*binaryCode*"0", "R"*binaryCode*"1")
	Sevaluator = (x) -> evaluator(x, "T"*binaryCode*"0", "S"*binaryCode*"1", "S"*binaryCode*"0", "Q"*binaryCode*"1")
	Revaluator = (x) -> evaluator(x, "R"*binaryCode*"0", "T"*binaryCode*"1", "Q"*binaryCode*"0", "R"*binaryCode*"1")
	Qevaluator = (x) -> evaluator(x, "R"*binaryCode*"0", "S"*binaryCode*"1", "Q"*binaryCode*"0", "Q"*binaryCode*"1")

	singValsT = bisectBetween(submatrixDict["T"*binaryCode*"0"], submatrixDict["T"*binaryCode*"1"], bk, Tevaluator, precB, limitB)
	singValsS = bisectBetween(submatrixDict["T"*binaryCode*"0"], submatrixDict["S"*binaryCode*"1"], bk, Sevaluator, precB, limitB)
	singValsR = bisectBetween(submatrixDict["R"*binaryCode*"0"], submatrixDict["T"*binaryCode*"1"], bk, Revaluator, precB, limitB)
	singValsQ = bisectBetween(submatrixDict["R"*binaryCode*"0"], submatrixDict["S"*binaryCode*"1"], bk, Qevaluator, precB, limitB)

	dict = Dict()
	push!(dict, "T"*binaryCode => singValsT)
	push!(dict, "S"*binaryCode => singValsS)
	push!(dict, "R"*binaryCode => singValsR)
	push!(dict, "Q"*binaryCode => singValsQ)

	return dict
end

""" Поиск собственных чисел родительской матрицы в интервалах между собственными числами дочерних матриц методом бисекций  
* eigvals1 - собственные числа первой матрицы;
* eigvals2 - собственные числа второй матрицы;
* bk - внедиагональный зануляемый при разделении родителькой матрицы элемент;
* evaluator - функция вычисляющая значение характерстического полинома;
* precB - критерий остановки метода бисекций, наибольшее отклонение характерстического полинома от нуля;
* limitB - критерий остановки метода бисекций, наибольшее количество делений пополам.""" 
function bisectBetween(eigvals1, eigvals2, bk, evaluator::Function, precB, limitB)
	eigVals = []
	unionSV = sort(vcat(eigvals1, eigvals2))
	len = length(unionSV)
	for i in 1:1:len - 1
		push!(eigVals, bisection(evaluator, unionSV[i], unionSV[i + 1], precB, limitB))
	end
	push!(eigVals, bisection(evaluator, unionSV[len], unionSV[len] + 2*bk, precB, limitB))

	return eigVals
end

function bisection(rootFunc::Function, a, b, prec, limit)
	rc = 0
    c = 0
    k = 1
    while(k < limit)
        k = k + 1
        c = (a + b) / 2
        ra = rootFunc(a)
        rb = rootFunc(b)
        rc = rootFunc(c)
        if (sign(ra) == sign(rb))
            return a
        end
        if (sign(ra) == sign(rc))
            a = c
        else
            b = c
        end
        if (abs(rc) < prec)
            return c
        end
    end
    return c
end

n = 100
println("n = ", n)

using PyPlot

setprecision(100)

godunovLists = lanczos(getGodunovMatrix(n)[1])
alphas = godunovLists[1]
bettas = godunovLists[2]

Af = toDense(alphas,bettas)
eigVal = sort(svdvals(Af), rev=true)

res = SingValsFinder(alphas, bettas, "", 1e-30, 1000)
a = sort(res["T"],rev=true)

a - eigVal


