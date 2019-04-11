@echo Generando objetos:
gfortran -c SCREEN3A.FOR
gfortran -c SCREEN3B.FOR
gfortran -c SCREEN3C.FOR
@echo Compilando SCREEN3:
gfortran -o screen3.exe SCREEN3A.o SCREEN3B.o SCREEN3C.o
@echo Compilación finalizada