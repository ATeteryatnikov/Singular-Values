""" функция возвращает косинус и синус - элементы матрицы Гивенса  """
function getGivensCS(a, b)
    r = sqrt(a^2+b^2)
    if (r==0)
        throw("Ошиба при вычислении чисел матрицы Гиаенса. Деление на ноль.")
    end
    c = a/r
    s = -b/r
    return [c,s]
end

function JacobiRotation(newAlphas, newBettas, numRepeat)
  
    if (length(newAlphas)==2)
        sv = svdvals(toDense(newAlphas,newBettas))
        return [sv, [0]]
    end

    N = length(newAlphas)  
    for numRep in 1:1:numRepeat
    # значения cos и sin для матрицы Гивенса (1,2).
        (cos0, sin0) = getGivensCS(newAlphas[1],newBettas[1])
    # после умножения трехдиагональной матрицы слева и справа на транспонированую матрицу Гивенса
    # полученная матрица теряет трехдиагональность из-за одного элемента выходящего за диагонали. 
    # Обозначим выступающий элемент 'P'.
    # Вычислим матричное умножение для симметричной матрицы по формулам полученным в mathcad:
        a1 = newAlphas[1]
        a2 = newAlphas[2]
        b1 = newBettas[1]
        b2 = newBettas[2]

        newAlphas[1] = a1*cos0^2 - 2*b1*cos0*sin0 + a2*sin0^2
        newBettas[1] = b1*cos0^2 - b1*sin0^2 + a1*cos0*sin0 - a2*cos0*sin0
        newAlphas[2] = a2*cos0^2 + 2*b1*cos0*sin0 + a1*sin0^2
        newBettas[2] = b2*cos0
        P = -b2*sin0

        for i in 1:1:length(newAlphas)-3
        # Для исключения выступающего элемента P вычисляем матрицу Гивенса (i+1,i+2)
            (cosN, sinN) = getGivensCS(newBettas[i], P)

            b1 = newBettas[i]
            b2 = newBettas[i+1]
            b3 = newBettas[i+2]
            a1 = newAlphas[i]
            a2 = newAlphas[i+1]
            a3 = newAlphas[i+2]
        # Умножаем слева и справа на матрицу Гивенса.
            newBettas[i] = b1*cosN - P*sinN
            newBettas[i+1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
            newBettas[i+2] = b3*cosN
            newAlphas[i+1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
            newAlphas[i+2] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
        # Получаем новый выступающий элемент P на позиции (i+2,i)       
            P = -b3*sinN
        end

    # Последняя итераци исключает P без появления нового

        (cosN, sinN) = getGivensCS(newBettas[N-2], P)

        b1 = newBettas[N-2]
        b2 = newBettas[N-1]
        a1 = newAlphas[N-2]
        a2 = newAlphas[N-1]
        a3 = newAlphas[N]

        newBettas[N-2] =  b1*cosN - P*sinN
        newBettas[N-1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
        newAlphas[N-1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
        newAlphas[N] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
    end

    return [newAlphas, newBettas]
end

function WilkinsonShift(lu,rd, b, i)
    delta = (rd - lu)/2
    shiftWilkinson = rd - (sign(delta)*b^2)/(abs(delta)+sqrt(delta^2 + b^2))

    return shiftWilkinson
end

function JacobiRotationWithShift(newAlphas, newBettas, numRepeat)
  
    if (length(newAlphas)==2)
        sv = svdvals(toDense(newAlphas,newBettas))
        return [sv, [0]]
    end

    N = length(newAlphas)  
    k=0
    firstElement = 1
    for numRep in 1:1:numRepeat
       
        #shiftWilkinson = WilkinsonShift(newAlphas[firstElement],newAlphas[firstElement+1],newBettas[firstElement+1], numRep)
        shiftWilkinson = WilkinsonShift(newAlphas[end-1],newAlphas[end],newBettas[end], numRep)
        #shiftWilkinson = GenericSVD.generic_svdvals!(toDense(alphas[end-1:end],bettas[end]))[2]
    # значения cos и sin для матрицы Гивенса (1,2).
        (cos0, sin0) = getGivensCS(newAlphas[firstElement] - shiftWilkinson, newBettas[firstElement])
    # после умножения трехдиагональной матрицы слева и справа на транспонированую матрицу Гивенса
    # полученная матрица теряет трехдиагональность из-за одного элемента выходящего за диагонали. 
    # Обозначим выступающий элемент 'P'.
    # Вычислим матричное умножение для симметричной матрицы по формулам полученным в mathcad:
       
        a1 = newAlphas[firstElement] - shiftWilkinson
        a2 = newAlphas[firstElement+1]- shiftWilkinson
        b1 = newBettas[firstElement]
        b2 = newBettas[firstElement+1]

        newAlphas[firstElement] = a1*cos0^2 - 2*b1*cos0*sin0 + a2*sin0^2 + shiftWilkinson
        newBettas[firstElement] = b1*cos0^2 - b1*sin0^2 + a1*cos0*sin0 - a2*cos0*sin0
        newAlphas[firstElement+1] = a2*cos0^2 + 2*b1*cos0*sin0 + a1*sin0^2 + shiftWilkinson
        newBettas[firstElement+1] = b2*cos0
        P = -b2*sin0
        
        for i in firstElement:1:length(newAlphas)-3
        # Для исключения выступающего элемента P вычисляем матрицу Гивенса (i+1,i+2)
           
            (cosN, sinN) = getGivensCS(newBettas[i], P)

            b1 = newBettas[i]
            b2 = newBettas[i+1]
            b3 = newBettas[i+2]
            a1 = newAlphas[i]
            a2 = newAlphas[i+1]
            a3 = newAlphas[i+2]
        # Умножаем слева и справа на матрицу Гивенса.
            newBettas[i] = b1*cosN - P*sinN
            newBettas[i+1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
            newBettas[i+2] = b3*cosN
            newAlphas[i+1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
            newAlphas[i+2] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2

        # Получаем новый выступающий элемент P на позиции (i+2,i)       
            P = -b3*sinN
        end

        # Последняя итераци исключает P без появления нового

        (cosN, sinN) = getGivensCS(newBettas[N-2], P)

        b1 = newBettas[N-2]
        b2 = newBettas[N-1]
        a1 = newAlphas[N-2]
        a2 = newAlphas[N-1]
        a3 = newAlphas[N]

        newBettas[N-2] =  b1*cosN - P*sinN
        newBettas[N-1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
        newAlphas[N-1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
        newAlphas[N] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
        newAlphas[N-2] = newAlphas[N-2]
    end
    return [newAlphas, newBettas, k]
end


function decomp(A,precision,limit=1000)
    A = copy(A)
    n = size(A)[1]
    Q = diagm(0=>ones(n))
    k=0 
    while(norm(diag(A,-1))>precision && k<limit)
        print(k,"\n")
        k=k+1
        delta = (A[n,n] - A[n-1,n-1])/2
        shiftWilkinson = A[n,n] - (sign(delta)*A[n,n-1]^2)/(abs(delta)+sqrt(delta^2 + A[n,n-1]^2))
        
        (Q1, R) = qr(A - shiftWilkinson*I)
        A = R*Q1 + shiftWilkinson*I
        #Q = Q1'*Q
    end
    #A = A + shiftWilkinson*I
    return A
end

function JacobiRotationRecursive(newAlphas, newBettas, numRepeat,lenBlock=4)

    k=0
    precDeflation=1e-100
    maxNorm=1e-20
    numRepeat=10

    unsv = []

    function deflationChecker(bettas, precDeflation)
        for i in 1:1:length(bettas)
            if (abs(bettas[i])<precDeflation)
                return i
            end
        end
        return -1
    end

    function recJacobi(newAlphas, newBettas, numRepeat, lenBlock)
        k=k+1
        if (k>10000)
            print("IterationEND")
            throw("Iterations END Limit")
            return
        end

        N = length(newAlphas) 
        if (N==1)
            unsv=union(unsv,newAlphas)
        elseif (N>1 && N<=lenBlock)
            # while(norm(betta)>maxNorm)
            #     alpha, betta = JacobiRotation(alpha,betta,numRepeat)
            # end
            #print("return by lenBlock with n: ",N,"\n")
            M = toDense(newAlphas, newBettas)
            a = svdvals(M)
            #decomp(M)
            unsv=union(unsv, a)
        elseif (N>lenBlock)
            alpha, betta = JacobiRotationWithShift(newAlphas,newBettas,numRepeat)
            if (norm(betta)<maxNorm)
                unsv=union(unsv,alpha)
                #print("return NULL\n")
                return
            end
            i = deflationChecker(newBettas,precDeflation)
            if (i>0)
                if (i==1)
                    #print("A ",i," ", N,"\n")
                    unsv=union(unsv,alpha[i])
                    recJacobi(alpha[i+1:end],betta[i+1:end],numRepeat,lenBlock)
                elseif (i>1 && i<N-1)
                    #print("B ",i," ", N,"\n")
                    recJacobi(alpha[1:i],betta[1:i-1],numRepeat,lenBlock)
                    recJacobi(alpha[i+1:end],betta[i+1:end],numRepeat,lenBlock)
                elseif (i==N-1)
                    #print("C ",i," ", N,"\n")
                    unsv=union(unsv,alpha[N])
                    recJacobi(alpha[1:i],betta[1:i-1],numRepeat,lenBlock)
                end
            else
                recJacobi(alpha,betta,numRepeat,lenBlock)
            end
        end
    end

    recJacobi(newAlphas, newBettas, numRepeat,lenBlock)
    
    return [unsv]
end



# include("src/Matrices.jl")
# include("Lanczos.jl")
# n = 200
# mantissa = 100
# numRepeated = 100
# maxBettasLimit = 1e-49
# lenBlock = 20

# setprecision(mantissa)
# #godunovLists = lanczos(rand(n,n))
# godunovLists = lanczos(getGodunovMatrix(n)[1])
# alphas = godunovLists[1]
# bettas = godunovLists[2]
# fullMatr = toDense(alphas, bettas)

# # JacobiRotation возвращает два массива (диагонали) и количество выполненых итераций [alphas, bettas, k]
# #JRRes, timerJR = @timed JacobiRotation(copy(alphas), copy(bettas), numRepeated)
# #JRRes, timerJR = @timed JacobiRotationModification(copy(alphas), copy(bettas), numRepeated, maxBettasLimit)
# originVals, denseTimer = @timed svdvals(fullMatr)
# jr, timerJR = @timed JacobiRotationWithShift(copy(alphas), copy(bettas), numRepeated)
# #jr2, timerJR = @timed JacobiRotation(copy(alphas), copy(bettas), numRepeated)
# vals, timerJRC = @timed JacobiRotationRecursive(copy(alphas), copy(bettas), numRepeated,lenBlock)
# vals = BigFloat.(vals[1])

# difference = sort(originVals,rev=true) - sort(abs.(jr[1]),rev=true)
# #difference2 = sort(originVals,rev=true) - sort(abs.(jr2[1]),rev=true)
# difference3 = sort(originVals,rev=true) - sort(abs.(vals),rev=true)
# maxBetta = maximum(jr[2])
# re = norm(difference)
# #re2 = norm(difference2)
# re3 = norm(difference3)