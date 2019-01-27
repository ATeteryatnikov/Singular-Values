
try
    type SparseNumber{}
        # номер строки
        row::Int
        # номера ненулевых элементов в строке
        notNullColumns::Array{Int}
        # номера строк исходного файла матрицы в которых хранится запись об элементах из notNullColumns
        rowsInSource::Array{Int}
    end
catch
    @showAndLog("SparseNumber уже подключен!")
end

""" Функция формирует массив элементы которого имеют тип SparseNumber"""
function createIndexStructure(Rows::Array{Int}, Columns::Array{Int})

    uniqueInd = sort!(unique(Rows))
    n = length(uniqueInd)
    ar = Array{SparseNumber}(n)
    for i in 1:1:n
        num = uniqueInd[i]
        rowNumberArray = findin(Rows, uniqueInd[i])
        notNullColumns = Columns[rowNumberArray]
        ar[i] = SparseNumber(num, notNullColumns, rowNumberArray)
    end
    return ar
end

""" jldFileName - название создаваемого файла jld;
TextSourceFile - массив значения BigFloat в текстовом виде;
IndStruct - Индексная структура SpareNumber"""
function writeInOneJLDGroup(jldFileName, TextSourceFile, IndStruct)

    BigFloatFunc = (x)->parse(BigFloat, x)
    file = jldopen(jldFileName, "w")

    for i in 1:1:length(IndStruct)
        stringValues = TextSourceFile[IndStruct[i].rowsInSource]
        bigFloatValues = map(BigFloatFunc, stringValues)
        file["BigFloatValues/$(IndStruct[i].row )"] = bigFloatValues
    end
    file.mmaparrays = true
    close(file)
end


""" Умножение разряженной матрицы на вектор """
function MulOnVec(jldFileWithValues, Struct, vector, tr=false)

    result = tr ? big.(zeros( width )) : big.(zeros( height ))
    for i in 1:1:size(Struct)[1]

        numRow = Struct[i].row
        column = vector[Struct[i].notNullColumns]
        row = read(jldFileWithValues["BigFloatValues/$(numRow)"])'

    	result[numRow] = ( row * column )[1]
    end
    return result
end
