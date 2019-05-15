@ECHO OFF
set dim=%1
set mantissa=%2
set prec=%3
set lenBlock=%4
set pathResults=%5
echo %pathResults%
set value=NormDiff
echo Agregate values for %value%*_dim_%dim%*_mantissa_%mantissa%*_bisectionPrec_%prec%*__lenBlock_%lenBlock%*
del "%pathResults%\graphicPoints.txt" 
for %%a in ("%pathResults%\%value%*_dim_%dim%*_mantissa_%mantissa%*_bisectionPrec_%prec%*_lenBlock_%numIterationlenBlock%*") do (
echo %%a
type "%%a"
echo. 
)>>"%pathResults%\graphicPoints.txt" 