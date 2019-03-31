"""
Функция вычисляет элементы матрицы вращений Гивенса.
* Принимает 2 числа:
    - a - элемент на главной диагонали;
    - b - зануляемый элемент.
* Возвращает массив из 2 элементов значения cos и sin
"""
function getGivensCS(a, b)
    r = sqrt(a^2 + b^2)
    if (r == 0)
        throw("Ошиба при вычислении чисел матрицы Гивенса. Деление на ноль.")
    end
    c = a / r
    s = -b / r
    return [c,s]
end

"""
Функция вычисляет сдвиг Уилкинсона для симметричной матрицы
по крайней правой нижней подматрицы размера 2x2.
* Принимает:
    - lu - левый верхний элемент подматрицы;
    - rd - правый нижний элемент подматрицы;
    - b - элемент вне главной диагонали подматрицы.
* Возвращает число - сдвиг Уилкинсона.
"""
function WilkinsonShift(lu,rd, b)
    delta = (rd - lu)/2
    shiftWilkinson = rd - (sign(delta)*b^2)/(abs(delta)+sqrt(delta^2 + b^2))

    return shiftWilkinson
end

""" 
Функция выполняет 'numIterations' итераций неявного QR-алгоритма
для симметричной трехдиагональной положительно определенной матрицы.
* Принимает:
    - alphas - главная диагональ;
    - bettas - диагональ соседняя с главной;
    - useShift - параметр использования сдвигов Уилкинсона;
    - numIterations - количество выполняемых итераций.
* Возвращает два массива - главную и соседнюю с главной диагонали.
"""
function ImplicitQRAlgorithm(alphas, bettas, useShift=false, numIterations=1)
    
    N = length(alphas)  

    for iteration in 1:1:numIterations
    # значения cos и sin для матрицы Гивенса (1,2).
        
        if (useShift) 
            shiftWilkinson = WilkinsonShift(alphas[end-1], alphas[end], bettas[end])
        else
            shiftWilkinson = 0
        end

        (cos0, sin0) = getGivensCS(alphas[1] - shiftWilkinson, bettas[1])
    # после умножения трехдиагональной матрицы слева и справа на транспонированую матрицу Гивенса
    # полученная матрица теряет трехдиагональность из-за одного элемента выходящего за диагонали. 
    # Обозначим выступающий элемент 'P'.
    # Вычислим матричное умножение для симметричной матрицы по формулам полученным в mathcad:
        
        a1 = alphas[1]
        a2 = alphas[2]
        b1 = bettas[1]
        b2 = bettas[2]
       
        alphas[1] = a1 * cos0^2 - 2 * b1 * cos0 * sin0 + a2 * sin0^2
        bettas[1] = b1 * cos0^2 - b1 * sin0^2 + a1 * cos0 * sin0 - a2 * cos0 * sin0
        alphas[2] = a2 * cos0^2 + 2 * b1 * cos0 * sin0 + a1 * sin0^2
        bettas[2] = b2 * cos0
        P = -b2 * sin0
        
        for i in 1:1:length(alphas) - 3
        # Для исключения выступающего элемента P вычисляем матрицу Гивенса (i+1,i+2)
            (cosN, sinN) = getGivensCS(bettas[i], P)

            b1 = bettas[i]
            b2 = bettas[i + 1]
            b3 = bettas[i + 2]
            a1 = alphas[i]
            a2 = alphas[i + 1]
            a3 = alphas[i + 2]
        # Умножаем слева и справа на матрицу Гивенса.
            bettas[i] = b1 * cosN - P * sinN
            bettas[i + 1] = b2 * cosN^2 - b2 * sinN^2 + a2 * cosN * sinN - a3 * cosN * sinN
            bettas[i + 2] = b3 * cosN
            alphas[i + 1] = a2 * cosN^2 - 2 * b2 * cosN * sinN + a3 * sinN^2
            alphas[i + 2] = a3 * cosN^2 + 2 * b2 * cosN * sinN + a2 * sinN^2
        # Получаем новый выступающий элемент P на позиции (i+2,i)       
            P = -b3 * sinN
        end
        
    # Последняя итераци исключает P без появления нового

        (cosN, sinN) = getGivensCS(bettas[N - 2], P)

        b1 = bettas[N - 2]
        b2 = bettas[N - 1]
        a1 = alphas[N - 2]
        a2 = alphas[N - 1]
        a3 = alphas[N]

        bettas[N - 2] =  b1 * cosN - P * sinN
        bettas[N - 1] = b2 * cosN^2 - b2 * sinN^2 + a2 * cosN * sinN - a3 * cosN * sinN
        alphas[N - 1] = a2 * cosN^2 - 2 * b2 * cosN * sinN + a3 * sinN^2
        alphas[N] = a3 * cosN^2 + 2 * b2 * cosN * sinN + a2 * sinN^2
    end

    return [alphas, bettas]
end

""" 
Функция выполняет 'numIterations' итераций неявного QR-алгоритма
для симметричной трехдиагональной положительно определенной матрицы.
* Принимает:
    - alphas - главная диагональ;
    - bettas - диагональ соседняя с главной;
    - prec - наибольшее отклонение по модулю нормы внедиагональных элементов от нуля;
    - useShift - параметр использования сдвигов Уилкинсона;
    - numQRIterations - количество выполняемых итераций между выполнением проверок
на наличие элемента из соседней с главной диагонали меньше чем prec; 
* Возвращает два массива - главную и соседнюю с главной диагонали.
"""
function RecursiveImplicitQR(alphas, bettas, prec, numQRIterations, useShift=false)

    # unsv - множество всех найденных сингулярных чисел
    unsv = []

    # Поиск индекса первого элемента из соседней с главной диагонали меньше prec
    function deflationChecker(bettas, prec)
        for i in 1:1:length(bettas)
            if (abs(bettas[i]) < prec)
                return i
            end
        end
        return -1
    end

    # Функция разделяет матрицу на две в случае, если в bettas присутствует элемент меньше prec
    function RecursiveImplicitQRDivider(alpha, betta, prec, numQRIterations, useShift=false)

        prevNormBetta = 0

        while(true)
            N = length(alpha) 
            if (N == 1)
                unsv = union(unsv, alpha)
                break
            elseif (N == 2)
            M = toDense(alpha, betta)
                unsv = union(unsv, svdvals(M))
                break
            else
                alpha, betta = ImplicitQRAlgorithm(alpha, betta, useShift, numQRIterations)
                normBetta = norm(betta)
                if (normBetta < prec || abs(prevNormBetta - normBetta) < prec)
                    unsv = union(unsv, alpha)
                    break
                end
                prevNormBetta = normBetta

                i = deflationChecker(betta, prec)
                if (i > 0)
                    if (i == 1)
                        unsv = union(unsv, alpha[i])
                        RecursiveImplicitQRDivider(
                            alpha[i + 1:end], betta[i + 1:end], prec, numQRIterations, useShift)
                        break
                    elseif (i > 1 && i < N - 1)
                        RecursiveImplicitQRDivider(
                            alpha[1:i], betta[1:i - 1], prec, numQRIterations, useShift)
                        RecursiveImplicitQRDivider(
                            alpha[i + 1:end], betta[i + 1:end], prec, numQRIterations, useShift)
                        break
                    elseif (i == N - 1)
                        unsv = union(unsv, alpha[N])
                        RecursiveImplicitQRDivider(
                            alpha[1:i], betta[1:i - 1], prec, numQRIterations, useShift)
                        break
                    end
                end
            end
        end
    end

    RecursiveImplicitQRDivider(
        alphas, bettas, prec, numQRIterations, useShift)
    
    return [unsv]
end