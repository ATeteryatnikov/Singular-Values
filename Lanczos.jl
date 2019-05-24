
try
    using GenericSVD
    using LinearAlgebra
catch ex
    print("Что-то пошло не так: ", ex)
end

""" 
Метод Ланцоша. Источник
LioydTtrefethen - Numerical Line 
page 277 (Modified for not Symmetric) 
Функция принимает:
	- A - прямоугольная вещественная матрица.
Функция возвращает:
	- alpha - главная диагональ симметричной трездиагональной матрицы;
	- betta - соседняя с главной диагональ.
"""
function lanczos(A)

    m = size(A)[1]
	n = size(A)[2]
    ### alpha - diagonal elements
   	alpha = big.(zeros(n))
   	### betta - subdiagonal elements
   	betta = big.(zeros(n-1))

    b = big.(ones(n))
    q = b/norm(b)

    q_prev = big.(zeros(n))
    a = big(0.0)
    b = big(0.0)

    for k in 1:1:n-1
        v = A'*(A*q)
        alpha[k] = q'*v

        v = v - b*q_prev - alpha[k]*q
        
        b = norm(v)
        
        betta[k] = b
        q_prev = q
        q = v/b
    end
    v = A'*(A*q)
    alpha[n] = q'*v

    return [alpha, betta]
end

function toDense(alpha, betta)

    len = length(alpha)
    A = fill(big(0e0), len, len)

    for i in 1:1:(len-1)
        A[i,i]=alpha[i]
        A[i+1,i] = A[i,i+1] = betta[i]
    end
    A[end,end] = alpha[end]

    return A
end
