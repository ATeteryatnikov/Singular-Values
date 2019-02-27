using DelimitedFiles

dir1 = "resultImplicitQR"
field = "DenseTime"
dim = "300"
mantissa = "*"
numIterations = "200"
prec="1e-80"

dir = pwd() * "\\" * dir1

filesList = readdir(dir)

filesList = filter(a->(findfirst(field,a)!=nothing),filesList)
filesList = filter(a->(findfirst(string(dim),a)!=nothing),filesList)
#filesList = filter(a->(findfirst(string(mantissa),a)!=nothing),filesList)
filesList = filter(a->(findfirst(string(numIterations),a)!=nothing),filesList)
filesList = filter(a->(findfirst(string(prec),a)!=nothing),filesList)

times = []
for fileName in filesList
	firstNumInFile = readdlm(dir * "\\" * fileName)[1]
	push!(times, firstNumInFile)
	print(dir * "\\" * fileName*"\n")
	print("time = ",firstNumInFile,"\n")
end