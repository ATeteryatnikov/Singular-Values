using DelimitedFiles
using PyPlot

availableColors = ["red","blue", "green", "brown", "purple"]

# параметры указывают, какие файлы попадают в выборку
dim = "100"
mantissa = "100"
rep = "*"

# Название графика в легенде
#label = string("dim=", dim, " mantissa=", mantissa, " rep=", rep)
label = string("dim=", dim, " mantissa=", mantissa)
color = rand(availableColors)
dirName = "resultJacobiRotationWithShiftTests"
createResultsBatPath = string(pwd(), "\\", "createResults.bat")
createResultsfilesPath = string(pwd(), "\\", dirName)
run(`cmd /c $createResultsBatPath $dim $mantissa $rep $createResultsfilesPath`)

# keyName определяет, какой параметр откладывается на ось абсцисс
keyName = "rep"
fileName = "graphicPoints.txt"
file = open(string(pwd(), "\\", dirName, "\\", fileName),"r")

lines = readlines(file)
n = length(lines)
close(file)
xMas = []
yMas = []

# Файл формируется чередованием пути к файлу и значением из файла
# в пути к файлу указаны параметры
for i in 1:1:n
	if (i%2!=0)
		ind1 = findlast(keyName, lines[i])[end]
		subStr = lines[i][ind1+2:end]
		ind2 = findfirst("_", subStr)
		str = subStr[1:ind2[1]-1]
		push!(xMas, parse(Int,str))
	else
		num = parse(BigFloat, lines[i])
		push!(yMas, num)
	end
end

indices = sortperm(xMas)
xMas = xMas[indices]
yMas = yMas[indices]

xlabel("Количество итераций")
ylabel("Норма разности сингулярных чисел")

plot(xMas,Float64.(yMas), "ro", label=label, color=color)
legend()
grid()