function norma(x)
  result = big(0.0)                       # Переменной result присваивается
                                          # тип BigFloat.
  x = x.^2         # Все элементы вектора возводятся в квадрат.
  result = sum(x)  # Все элементы вектора суммируются.
  return result    # функция возвращает переменную result.
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
