
"""макрос выводит информацию на дисплей и в лог-файл. предполагает, что logFileName и path определены до использования макроса"""
function showAndLog(str)
    logfile = open(string(path, logFileName), "a+")
    str = string(eval(str), "\n")
    write(logfile, str)
    showall(str)
    close(logfile)
end
