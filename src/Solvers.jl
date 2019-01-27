function norma(x)
  x = x.^2
  result = sum(x)
  return result
end

function doMCG(A,B,x, norm_residual, limit)
  p = A'*(B-A*x)
  s=p
  y=norma(s)
  k = 0
  while(y> norm_residual && k<limit)
    q = A*p
    a = y/norma(q)
    x = x + a*p
    s = s - a*A'*q
    y1 = norma(s)
    b = y1/y
    p = s + b*p
    k=k+1
    y=y1
  end
  return(x)
end


""" QR Decomposition. function get (A) - square matrix """
function doQRDecomposition(A)

    R = copy(A)
    last = size(A)[1]
    Q = eye(last, last)
    for i in 1:1:(last-1)

        column_i = R[i:end, i]

        alpha = norm(column_i)

        delta_i = -column_i
        delta_i[1]=delta_i[1] + alpha

        n = delta_i / norm(delta_i)

        for k in (i+1):1:last
            R[i:end,k] = R[i:end,k] - 2*n*(R[i:end,k]'*n)
        end
        R[i,i] = alpha
        R[i+1:last,i] = zeros(last-i)

        for k in 1:1:last
            Q[i:end,k] = Q[i:end,k] - 2*n*(Q[i:end,k]'*n)
        end
    end

    return [Q',R]
end

""" Swap 2 columns in Matrix """
function TranspositionColumns(A, i, j)
    temp = A[:,i]
    A[:,i] = A[:,j]
    A[:,j] = temp
end

""" Swap 2 rows in Matrix """
function TranspositionRows(A, i, j)
    temp = A[i,:]
    A[i,:] = A[j,:]
    A[j,:] = temp
end

""" LU - Decomposition. ( P*A*Q = LU ) """
function DoLUDecomposition(A)
    last = size(A)[1]
    CompactMatrix = copy(A)
    P = eye(last, last)
    Q = eye(last, last)

    for diagNum in 1:1:last

        num = findmax(abs(CompactMatrix[diagNum:end, diagNum:end]))[2]

        if (num == 0)
            print("\n Ошибка \n")
        elseif (mod(num, last-diagNum+1)==0)
            numCol = diagNum + div(num, last-diagNum+1)-1
            numlast = last
        elseif (div(num, last-diagNum+1)==0)
            numCol = diagNum + 1
            numlast = diagNum + mod(num, last-diagNum+1)
        else
            numCol = diagNum + div(num, last-diagNum+1)
            numlast = diagNum + mod(num, last-diagNum+1) - 1
        end

        TranspositionColumns(CompactMatrix, numCol, diagNum)
        TranspositionColumns(Q, numCol, diagNum)
        TranspositionRows(CompactMatrix, numlast, diagNum)
        TranspositionRows(P, numlast, diagNum)

        for i in (diagNum+1):1:last
            for j = (diagNum+1):1:last
                CompactMatrix[i,j] = CompactMatrix[i,j] - CompactMatrix[diagNum, j] * CompactMatrix[i,diagNum] / CompactMatrix[diagNum,diagNum]
            end
        end
        for j = (diagNum+1):1:last
            CompactMatrix[diagNum,j] = CompactMatrix[diagNum,j]/CompactMatrix[diagNum,diagNum]
        end
    end
    L = copy(LowerTriangular(CompactMatrix))
    U = copy(UpperTriangular(CompactMatrix))
    for i in 1:1:5 U[i,i]=1 end

    return [P,L,U,Q]    #  P'*L*U*Q'
end
