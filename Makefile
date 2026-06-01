CXX    := g++
AR     := ar
LDFLAGS := -lopenblas -lm -lpthread

BUILD ?= release
ifeq ($(BUILD),debug)
  CFLAGS := -std=c++17 -O0 -g -Wall -Wextra -Iinclude -MMD -MP
else
  CFLAGS := -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Iinclude -MMD -MP
endif

BLAS ?= 1
ifeq ($(BLAS),1)
  CFLAGS += -DUSE_BLAS
endif

VARIANTE ?= $(if $(filter 1,$(BLAS)),blas,escalar)
CFLAGS   += -DVARIANTE='"$(VARIANTE)"'

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
	rm -rf build $(LIB) $(PLOT_CSV)

# --- Tests ---
# Lista todos os .cpp em tests/ e define os binários correspondentes em build/tests/
TEST_SRCS := $(filter-out tests/benchmark.cpp, $(wildcard tests/*.cpp))
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

COMPARE_CSV := plot/comparar.csv
COMPARE_PNG := plot/comparar.png

# Compila e roda as variantes blas e escalar, mescla os CSVs e gera o gráfico
comparar:
	$(MAKE) clean && $(MAKE) $(BENCHMARK_BIN) BLAS=1 VARIANTE=blas
	./$(BENCHMARK_BIN) > $(COMPARE_CSV)
	$(MAKE) clean && $(MAKE) $(BENCHMARK_BIN) BLAS=0 VARIANTE=escalar
	./$(BENCHMARK_BIN) | tail -n +2 >> $(COMPARE_CSV)
	$(PYTHON) plot/plot.py $(COMPARE_CSV) $(COMPARE_PNG)

PYTHON    := plot/.venv/bin/python
PLOT_CSV  := plot/data.csv

plot/%.png: $(BENCHMARK_BIN) plot/plot.py
	./$(BENCHMARK_BIN) > $(PLOT_CSV)
	$(PYTHON) plot/plot.py $(PLOT_CSV) $@

# --- Instalação ---
INSTALL_LIB_DIR  := /usr/local/lib
INSTALL_INC_DIR  := /usr/local/include/linear_systems

install: $(LIB)
	mkdir -p $(INSTALL_LIB_DIR)
	mkdir -p $(INSTALL_INC_DIR)
	cp $(LIB) $(INSTALL_LIB_DIR)/
	cp -r include/* $(INSTALL_INC_DIR)/

uninstall:
	rm -f $(INSTALL_LIB_DIR)/$(LIB)
	rm -rf $(INSTALL_INC_DIR)

.PHONY: all clean test benchmark install uninstall
