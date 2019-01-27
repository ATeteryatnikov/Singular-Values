function norma(x)
  result = big(0.0)
  x = x.^2
  result = sum(x)
  return result
end

""" A - открытый на чтение jld файл с BigFloat значениями для матрицы
A_tr - открытый на чтение jld файл с BigFloat значениями для траспонированной матрицы
B - Правая часть СЛАУ """
function doMCG(A, A_Struct,
               A_tr, A_tr_Struct,
               B, x,  limit=1, norm_residual=0)

  Ax = MulOnVec(A, A_Struct, x)
  p = MulOnVec(A_tr, A_tr_Struct, (B-Ax), true)
  s=p
  y=norma(s)
  Ax = 0
  k = 0
  while(y> norm_residual && k<limit)
    q = MulOnVec(A, A_Struct, p)
    a = y/norma(q)
    x = x + a*p
    s = s - a * MulOnVec(A_tr, A_tr_Struct, q, true)
    y1 = norma(s)
    b = y1/y
    p = s + b*p
    k=k+1
    y=y1
  end
  return(x)
end
