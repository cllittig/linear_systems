#Makefile da biblioteca

#INICIO DAS VARIAVEIS DE AMBIENTE
#Definindo variaveis de ambiente .hpp
#=========================================================================

#Variaveis include/algebra_linear

MATRIZHPP = include/algebra_linear/matriz.hpp
VECTORHPP = include/algebra_linear/vetor.hpp
SISTEMA_LINEARHPP = include/algebra_linear/sistema_linear.hpp
#=========================================================================

#variaveis include/metodos

GAUSS_JORDANHPP = include/metodos/gauss_jordan.hpp
#=========================================================================

#variaveis include/utils

READ_DATAHPP = include/utils/read_data.hpp
#=========================================================================

#Definindo variaveis .cpp
#=========================================================================

#variavel main.cpp

MAIN = src/main.cpp

#variaveis src/algebra_linear

MATRIZCPP = src/algebra_linear/matriz.cpp
VECTORCPP = src/algebra_linear/vetor.cpp
SISTEMA_LINEARCPP = src/alebra_linear/sistema_linear.cpp
#=========================================================================

#variaveis src/metodos

GAUSS_JORDANCPP = src/metodos/gauss_jordan.cpp
#=========================================================================

READ_DATACPP = src/utils/read_data.cpp
#========================================================================

#FIM DAS VARIAVEIS DE AMBIENTE

#INICIO DAS REGRAS DE COMPILAÇÃO

#linka o código-objeto
main: main.o matriz.o gauss_jordan.o
	g++ -o main main.o matriz.o gauss_jordan.o

#compila o main.cpp
main.o: $(MAIN) $(MATRIZHPP)
	g++ -c -Iinclude -o main.o $(MAIN)

#compila matriz.cpp
matriz.o: $(MATRIZCPP) $(MATRIZHPP)
	g++ -c -Iinclude -o matriz.o $(MATRIZCPP)	

#compila gauss_jordan.cpp
gauss_jordan.o: $(GAUSS_JORDANCPP) $(GAUSS_JORDANHPP) $(MATRIZHPP)
	g++ -c -Iinclude -o gauss_jordan.o $(GAUSS_JORDANCPP) 

clean:
	rm -f main *.o
