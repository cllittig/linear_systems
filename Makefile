CXX    := g++
AR     := ar
# -Iinclude faz com que o compilador busque dentro da pasta include
CFLAGS := -std=c++17 -O3 -Wall -Wextra -Iinclude -MMD -MP
LDFLAGS := -lopenblas -lm -lpthread

# Busca todos os arquivos .cpp em src e suas subpastas
SOURCE := $(shell find src -name '*.cpp')
# Define onde os objetos .o ficarão
OBJS := $(SOURCE:src/%.cpp=build/%.o)
DEPS := $(OBJS:.o=.d)

LIB := libls.a

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

# Regra de compilação
build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) -c $< -o $@

# Limpeza
clean:
	rm -rf build $(LIB) $(PLOT_CSV) $(PLOT_PNG)

# --- Tests ---
# Lista todos os .cpp em tests/ e define os binários correspondentes em build/tests/
TEST_SRCS := $(wildcard tests/*.cpp)
TEST_BINS := $(TEST_SRCS:tests/%.cpp=build/tests/%)

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

# Regra genérica: compila tests/foo.cpp linkando contra a biblioteca
build/tests/%: tests/%.cpp $(LIB)
	@mkdir -p build/tests
	$(CXX) $(CFLAGS) $< -L. -lls $(LDFLAGS) -o $@

-include $(DEPS)

# --- Benchmark ---
BENCHMARK_BIN := build/benchmark

benchmark: $(BENCHMARK_BIN)
	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): tests/benchmark.cpp $(LIB)
	@mkdir -p build
	$(CXX) $(CFLAGS) $< -L. -lls $(LDFLAGS) -o $@

PYTHON    := plot/.venv/bin/python
PLOT_CSV  := plot/data.csv
PLOT_PNG  := plot/data.png

$(PLOT_PNG): $(BENCHMARK_BIN) plot/plot.py
	./$(BENCHMARK_BIN) > $(PLOT_CSV)
	$(PYTHON) plot/plot.py $(PLOT_CSV)

# --- Installation ---
INSTALL_LIB_DIR  := /usr/local/lib
INSTALL_INC_DIR  := /usr/local/include/linear_systems

install: $(LIB)
	@echo "Instalando biblioteca..."
	mkdir -p $(INSTALL_LIB_DIR)
	mkdir -p $(INSTALL_INC_DIR)
	cp $(LIB) $(INSTALL_LIB_DIR)/
	cp -r include/* $(INSTALL_INC_DIR)/
	@echo "✓ Instalado em $(INSTALL_LIB_DIR) e $(INSTALL_INC_DIR)"

uninstall:
	@echo "Removendo biblioteca..."
	rm -f $(INSTALL_LIB_DIR)/$(LIB)
	rm -rf $(INSTALL_INC_DIR)
	@echo "✓ Removido"

.PHONY: all clean test benchmark install uninstall
