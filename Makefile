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

# --- Tests ---
# Lista todos os .cpp em tests/ e define os binários correspondentes em build/tests/
TEST_SRCS := $(wildcard tests/*.cpp)
TEST_BINS := $(TEST_SRCS:tests/%.cpp=build/tests/%)

# Objetos da biblioteca sem o main do projeto — cada teste tem seu próprio main(),
# então linkar build/main.o causaria símbolo duplicado.
SRC_OBJS_NO_MAIN := $(filter-out build/main.o, $(OBJS))

# Compila e executa todos os testes. Retorna erro se algum falhar.
test: $(TEST_BINS)
	@echo ""; failed=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t ---"; $$t; [ $$? -ne 0 ] && failed=$$((failed+1)); \
	done; \
	echo ""; [ $$failed -eq 0 ] \
		&& echo "Todos os testes passaram." \
		|| echo "$$failed teste(s) falharam."; \
	exit $$failed

# Regra genérica: compila tests/foo.cpp linkando com os objetos da biblioteca
build/tests/%: tests/%.cpp $(SRC_OBJS_NO_MAIN)
	@mkdir -p build/tests
	$(CXX) $(CFLAGS) $< $(SRC_OBJS_NO_MAIN) $(LDFLAGS) -o $@

-include $(DEPS)

.PHONY: all clean run test
