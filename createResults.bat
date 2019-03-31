@ECHO OFF
set dim=%1
set mantissa=%2
set prec=%3
set numIteration=%4
set pathResults=%5
echo %pathResults%
set value=NormDiff
echo Agregate values for %value%*_dim_%dim%*_mantissa_%mantissa%*_prec_%prec%*__numIteration_%numIteration%*
del "%pathResults%\graphicPoints.txt" 
for %%a in ("%pathResults%\%value%*_dim_%dim%*_mantissa_%mantissa%*_prec_%prec%*_numIteration_%numIteration%*") do (
echo %%a
type "%%a"
echo. 
)>>"%pathResults%\graphicPoints.txt" 