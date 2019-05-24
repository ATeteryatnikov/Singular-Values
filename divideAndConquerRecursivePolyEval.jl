# Алгоритм divide-and-conquer, вычисление значения полинома производится рекурсивно для каждой точки
# "V. Rokhlin, A fast divide-and-conquer algorithm 
# for computing the spectra of real symmetric tridiagonal matrices"
include("src/LogGenerator.jl")
include("lanczos.jl")
include("src/Matrices.jl")     
include("implicitQR.jl")


"""
 Вычисление значения характеристического полинома в точке x по теореме 4.6 из статьи "V. Rokhlin
A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices"
Функция возвращает значение характеристического полинома в точке.
Функция принимает на вход:
* x - точка в которой вычисляется значение характеристического полинома;
* alphas - главная диагональ трехдиагональной матрицы;
* bettas - соседняя с главной диагональ трехдиагональной матрицы;
* lenBlock - максимальная размерность матриц в рекурсивном разделении матриц 
при которой расчет значения характеристического полинома происходит через 'det'. 
"""
function PolyRoot(x, alpha, betta, lenBlock, evalPsi=false)

    len = length(alpha)
    if (len <= lenBlock && len > 0)
        M = toDense(alpha, betta)
        return det(M - x * I)
    elseif (len == 0)
        return 1
    else
        k = div(len, 2)
        X1 = X2 = X3 = X4 = big(1.0)

        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k + 1] =  alpha[k + 1] - abs(betta[k])

        if (k > 0)
            X1 = PolyRoot(x, alpha[1:k], betta[1:k - 1], lenBlock)
        end
        if (len - k > 0)
            X2 = PolyRoot(x, alpha[k + 1:len], betta[k + 1:len - 1], lenBlock)
        end
        if (k - 1 > 0)
            X3 = PolyRoot(x, alpha[1:k - 1], betta[1:k - 2], lenBlock)
        end
		if (len - k - 1 > 0)
            X4 = PolyRoot(x, alpha[k + 2:len], betta[k + 2:len - 1], lenBlock)
        end

        alpha[k] = alpha[k] + abs(betta[k])
        alpha[k + 1] =  alpha[k + 1] + abs(betta[k])

        X0 = X1 * X2 + abs(betta[k]) * X1 * X4 + abs(betta[k]) * X2 * X3
        
        fi0 = -log(abs(betta[k])) - log(abs(X4/X2 + X3/X1))
        return evalPsi ? fi0 : X0 

        #return X0
    end
end

"""
Вычисление собственных чисел плотной матрицы QR-алгоритмом
"""
function denseQRAlgorithm(A, prec)
    A = copy(A)
    Q = Diagonal(ones(size(A)[1]))
    r=1
    while (r>prec)
        (Q1, R) = qr(A)
        A = R * Q1
        #Q = Q1'*Q
        r = norm(diag(A,-1))
    end
return diag(A)
    #return [Q, A]
end

"""
 Вычисление корня характеристического полинома методом половинного 
деления на интервале.
Вычисление значений производится функцией PolyRoot.
Функция возвращает корень характеристического полинома матрицы.
Функция принимает на вход:
* a, b - точки начала и конца интервала в котором происходит поиск корня;
* alphas - главная диагональ трехдиагональной матрицы;
* bettas - соседняя с главной диагональ трехдиагональной матрицы;
* lenBlock - (для метода PolyRoot) Максимальная размерность матриц в рекурсивном 
разделении матриц при которой расчет значения характеристического 
полинома происходит через 'det';
* prec - при нахождении значения характеристического полинома меньше чем prec, 
считаем, что корень найден; 
* limit - ограничение количества итераций. 
"""
function bisection(a, b, alphas, bettas, lenBlock = 4, prec = 1e-15, limit = 1000)
    root = x->PolyRoot(x, (alphas), (bettas), lenBlock)
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
            ra = root(a)
            a_changed = false
        end
        if (b_changed)
            rb = root(b)
            b_changed = false
        end
        
        if (sign(ra) == sign(rb))
            return a
        end
        rc = root(c)

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

"""
 Алгоритм "divide-and-conquer" 
(описывается в пукте 4 статьи "V. Rokhlin A fast divide-and-conquer algorithm for computing
 the spectra of real symmetric tridiagonal matrices")
Рекурсивный алгоритм вычисления собственных чисел.
Функция возвращает массив собственных чисел.
Функция принимает на вход:
* alphas - главная диагональ трехдиагональной матрицы;
* bettas - соседняя с главной диагональ трехдиагональной матрицы;
* lenBlock - (для метода PolyRoot) Максимальная размерность матриц в рекурсивном разделении 
матриц при которой расчет значения характеристического полинома происходит через 'det', а 
расчет собственных чисел через 'denseQRAlgorithm';
* prec - при нахождении значения характеристического полинома меньше чем prec, считаем, 
что корень найден; 
* step - Отладочная информация, позволяет определить, какая часть матрицы используется. 
"""
function svdValsFinder(alpha, betta, lenBlock, prec = 1e-15, step = "")

    len = length(alpha)
    #println("step: ",step)
    if (len <= lenBlock)
        return RecursiveImplicitQR(copy(alpha), copy(betta), prec, 200)
        #matr = toDense(alpha, betta)
        #return denseQRAlgorithm(matr, prec)
    else
        k = div(len, 2)

        alphas = copy(alpha)
        bettas = copy(betta)
        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k + 1] = alpha[k + 1] - abs(betta[k])

        if (k - 1 >= 0)
            singValsT0 = svdValsFinder((alpha[1:k]), (betta[1:k - 1]), lenBlock, prec, step * "0")
        end
        if (len - (k + 1) >= 0)
            singValsT1 = svdValsFinder((alpha[k + 1:end]), (betta[k + 1:end]), lenBlock, prec, step * "1")
        end

        alpha[k] = alpha[k] + abs(betta[k])
        alpha[k + 1] = alpha[k + 1] + abs(betta[k])

        unionSingVals = sort(vcat(singValsT0, singValsT1))
        
        singVals = []
        for i in 1:1:len - 1
            #println("start bisection ", i, " of ", len)
            push!(singVals, bisection(unionSingVals[i], unionSingVals[i + 1], alpha, betta, lenBlock, prec))
        end
        push!(singVals, bisection(unionSingVals[len], unionSingVals[len] + 2 * betta[k], alpha, betta, lenBlock, prec))

        return singVals
    end
end



