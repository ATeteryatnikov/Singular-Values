using DelimitedFiles

dir1 = "resultJacobiRotationWithoutShiftTests"
field = "JRTime"
dim = 300
mantissa = 100

dir = pwd() * "\\" * dir1

filesList = readdir(dir)

filesList = filter(a->(findfirst(field,a)!=nothing),filesList)
filesList = filter(a->(findfirst(string(dim),a)!=nothing),filesList)
filesList = filter(a->(findfirst(string(mantissa),a)!=nothing),filesList)

times = []
for fileName in filesList
	print(dir * "\\" * fileName*"\n")
	push!(times, readdlm(dir * "\\" * fileName)[1])
end