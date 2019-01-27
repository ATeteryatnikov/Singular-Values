@ECHO OFF
set dim=%1
set mantissa=%2
set rep=%3
set pathResults=%4
echo %pathResults%
set value=NormDiff
echo Agregate values for %value%*_dim_%dim%*_mantissa_%mantissa%*_rep_%rep%*
del "%pathResults%\graphicPoints.txt" 
for %%a in ("%pathResults%\%value%*_dim_%dim%*_mantissa_%mantissa%*_rep_%rep%*") do (
echo %%a
type "%%a"
echo. 
)>>"%pathResults%\graphicPoints.txt" 