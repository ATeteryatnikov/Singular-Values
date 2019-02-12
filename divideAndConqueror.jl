include("src/LogGenerator.jl")
include("lanczos.jl")
include("src/Matrices.jl")     

""" Вычисление полинома в точке x по схеме  4.6 из статьи V. Rokhlin
A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices """
function PolyRoot(x, alpha, betta, lenBlock, evalPsi=false)

    getDet = (alpha, betta, x) -> det(toDense(alpha, betta)-x*I)

    len = length(alpha)
    if (len <= lenBlock && len>0)
        M = toDense(alpha, betta)
        return det(M-x*I)
    elseif (len==0)
        return 1
    else
        k = div(len, 2)
        X1 = X2 = X3 = X4 = big(1.0)

        a1 = copy(alpha[k])
        a2 = copy(alpha[k+1])
        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k+1] =  alpha[k+1] - abs(betta[k])

        if (k > 0)
            X1 = PolyRoot(x,alpha[1:k], betta[1:k-1], lenBlock)
        end
        if (len-k > 0)
            X2 = PolyRoot(x, alpha[k+1:len], betta[k+1:len-1], lenBlock)
        end
        if (k-1 > 0)
            X3 = PolyRoot(x, alpha[1:k-1], betta[1:k-2], lenBlock)
        end
        if (len-k-1 > 0)
            X4 = PolyRoot(x, alpha[k+2:len], betta[k+2:len-1], lenBlock)
        end

        alpha[k] = alpha[k] + abs(betta[k])
        alpha[k+1] =  alpha[k+1] + abs(betta[k])

        X0 = X1*X2 + abs(betta[k])*X1*X4 + abs(betta[k])*X2*X3
        fi0 = -log(abs(betta[k])) - log(abs(X4/X2 + X3/X1))
        return evalPsi ? fi0 : X0 
        return X0
    end
end

"Вычисление сингулярных чисел плотной матрицы QR-алгоритмом"
function SVDQR(A)
    A = copy(A)
    Q = Diagonal(ones(size(A)[1]))
    #Q1 = Diagonal(ones(size(A)[1]))
    for i in 1:1:100
        (Q1, R) = qr(A)
        A = R*Q1
        #Q = Q1'*Q
    end
    return diag(A)
    #return [Q, A]
end

function bisection(a,b, alphas, bettas, lenBlock=4, prec=1e-15)
    root = x->PolyRoot(x, (alphas), copy(bettas),lenBlock)
    stepLimit=1000
    c=0
    k=1
    while(k<stepLimit)
        k=k+1
        c = (a+b)/2
        ra = root(a)
        rb = root(b)
        rc = root(c)
        if (sign(ra)==sign(rc))
            a = c
        else
            b = c
        end
        if (abs(rc)<prec)
            return c
        end
    end
    return c
end

""" Рекурсивный поиск сингулярных чисел между сингулярными числами дочерных матриц (Т0 и Т1) """
function svdValsFinder(alpha, betta, lenBlock, prec=1e-15, step = 0)

    len = length(alpha)
    step = step+1
    if (len <= lenBlock)
        matr = toDense(alpha, betta)
        return SVDQR(matr)
    else
        k = div(len, 2)

        alphas = copy(alpha)
        bettas = copy(betta)
        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k+1] = alpha[k+1] - abs(betta[k])

        if (k-1 >= 0)
            singValsT0 = svdValsFinder(copy(alpha[1:k]), copy(betta[1:k-1]), lenBlock, prec,step)
        end
        if (len-(k+1) >= 0)
            singValsT1 = svdValsFinder(copy(alpha[k+1:end]), copy(betta[k+1:end]), lenBlock, prec, step)
        end

        alpha[k] = alpha[k] + abs(betta[k])
        alpha[k+1] = alpha[k+1] + abs(betta[k])

        unionSingVals = sort(vcat(singValsT0, singValsT1))
        
        singVals = []
        for i in 1:1:len-1
            push!(singVals, bisection(unionSingVals[i], unionSingVals[i+1], alphas, bettas, lenBlock,prec))
        end
        push!(singVals, bisection(unionSingVals[len], unionSingVals[len]+2*betta[k], alphas, bettas, lenBlock, prec))

        return singVals
    end
end



