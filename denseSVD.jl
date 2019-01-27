include("Lanczos.jl")
include("JacobiRotation.jl")
include("src/Matrices.jl")
using GenericSVD
using DelimitedFiles

mantissa = 100
dim = 25

setprecision(mantissa)
godunovM = getGodunovMatrix(dim)[1]
# alphas = godunovLists[1]
# bettas = godunovLists[2]
# fullMatr = toDense(alphas, bettas)
singVals = svdvals(godunovM)


function decomp(A)
    A = copy(A)
    n = size(A)[1]
    Q = diagm(0=>ones(n))
    ss = 0

    k=n
    for i in 2:1:1000
    	
    	if (abs(A[k-1,k])>1e-20)
	    	delta = (A[k-1,k-1] - A[k,k])/2
			shiftWilkinson = A[k,k] - (sign(delta)*A[k-1,k]^2)/(abs(delta)+sqrt(delta^2 + A[k-1,k]^2))
			#shiftWilkinson = A[n,n] + delta - (sign(delta)*sqrt(delta^2 + A[n,n-1]^2))
			#shiftWilkinson = A[end,end]
			
			println("shift = ",shiftWilkinson, " , k = ", k ,  "\n")
		else
			shiftWilkinson = 0
			if (k>2)
				k = k - 1
			else
				k=n
			end
		end
		#shiftWilkinson=0
        (Q1, R) = qr(A - shiftWilkinson*I)
        A = R*Q1 + shiftWilkinson*I
        Q = Q1'*Q
        if (i%100 == 99)

        	println("shift = ",shiftWilkinson, "\n")
        	re = sqrt.(sort(abs.(diag(A)),rev=true)) - singVals
        	plot(i,Float64(norm(re)),"ro")
        end
    end
    #A = A + shiftWilkinson*I
    return [Q, A]
end



A = decomp(godunovM'*godunovM)[2]
re = sqrt.(sort(abs.(diag(A)),rev=true)) - singVals
nre = norm(re)