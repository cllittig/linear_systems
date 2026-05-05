CXX := g++
# -Iinclude faz com que o compilador busque dentro da pasta include
CFLAGS := -std=c++17 -O3 -Wall -Wextra -Iinclude -MMD -MP
# LDFLAGS corrigido e pronto para o linker
LDFLAGS := -lopenblas -lm -lpthread

# Busca todos os arquivos .cpp em src e suas subpastas
SOURCE := $(shell find src -name '*.cpp')
# Define onde os objetos .o ficarão
OBJS := $(SOURCE:src/%.cpp=build/%.o)
DEPS := $(OBJS:.o=.d)

# Nome do executável final
TARGET := programa_tcc

# Alvo principal: compila e gera o executável
all: $(TARGET)

# Linkagem: LDFLAGS DEVE ir ao final da linha
$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

# Regra de Compilação: Aqui é onde o CFLAGS (com o -Iinclude) entra
build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) -c $< -o $@

# Atalho para rodar
run: all
	./$(TARGET)

# Limpeza
clean:
	rm -rf build $(TARGET)

-include $(DEPS)

.PHONY: all clean run
