# COMPILAR PARA LINUX
FFLAGS = 
OBJS = SCREEN3A.o SCREEN3B.o SCREEN3C.o

all: start screen3 end

start:
	@echo Generando objetos:

screen3: $(OBJS)
	@echo Compilando SCREEN3:
	gfortran -o $@ $?

%.o: %.f
	gfortran -c $< -o $@

end:
	@echo Compilación finalizada.

clean:
	rm *.o
