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

LAPACK ?= 0
ifeq ($(LAPACK),1)
  CFLAGS  += -DUSE_LAPACK -I/usr/include/openblas
endif

VARIANTE ?= $(if $(filter 1,$(LAPACK)),lapack,$(if $(filter 1,$(BLAS)),blas,escalar))
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
# Google Test Configuration
GTEST_LIBS := -lgtest -lgtest_main
GTEST_CPPFLAGS := $(shell pkg-config --cflags gtest 2>/dev/null)

# Octave Configuration (for comparative tests)
# Octave >= 6 removeu os flags octave-config --cppflags/--libs/--octlibdir.
# A interface estável entre versões e o `mkoctfile -p` para flags de compilação
# e `octave-config -p VAR` para diretórios.
OCTAVE_OCTLIBDIR := $(shell octave-config -p OCTLIBDIR 2>/dev/null)
OCTAVE_INCFLAGS  := $(shell mkoctfile -p INCFLAGS 2>/dev/null)
OCTAVE_LIBS      := $(shell mkoctfile -p OCTAVE_LIBS 2>/dev/null)

# Test compilation flags
TEST_CPPFLAGS := $(CFLAGS) $(GTEST_CPPFLAGS) $(OCTAVE_INCFLAGS)
TEST_LDFLAGS := $(LDFLAGS) $(GTEST_LIBS)
ifneq ($(OCTAVE_OCTLIBDIR),)
  # -rpath garante que o binário encontre a liboctave.so em tempo de execução
  # (ela fica fora do caminho padrão do linker dinâmico).
  TEST_LDFLAGS += -L$(OCTAVE_OCTLIBDIR) -Wl,-rpath,$(OCTAVE_OCTLIBDIR) $(OCTAVE_LIBS)
endif

# Lista todos os .cpp em tests/ e define os binários correspondentes em build/tests/
TEST_SRCS := $(filter-out tests/benchmark.cpp tests/benchmark_runner.cpp tests/test_octave_comparison.cpp tests/octave_ref.cpp, $(wildcard tests/*.cpp))
TEST_BINS := $(TEST_SRCS:tests/%.cpp=build/tests/%)

# Octave comparison test (only if Octave is available)
TEST_OCTAVE_SRC := tests/test_octave_comparison.cpp
TEST_OCTAVE_BIN := build/tests/test_octave_comparison
ifneq ($(OCTAVE_OCTLIBDIR),)
  TEST_BINS_WITH_OCTAVE := $(TEST_BINS) $(TEST_OCTAVE_BIN)
else
  TEST_BINS_WITH_OCTAVE := $(TEST_BINS)
endif

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

# Compila e executa testes incluindo comparação com Octave
test-octave: $(TEST_BINS_WITH_OCTAVE)
	@echo ""; failed=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t ---"; $$t; [ $$? -ne 0 ] && failed=$$((failed+1)); \
	done; \
	if [ -f "$(TEST_OCTAVE_BIN)" ]; then \
		echo ""; echo "=== Executando testes comparativos com Octave ==="; \
		$(TEST_OCTAVE_BIN); [ $$? -ne 0 ] && failed=$$((failed+1)); \
	else \
		echo ""; echo "WARNING: Octave comparison tests not available (Octave not found)"; \
	fi; \
	echo ""; [ $$failed -eq 0 ] \
		&& echo "Todos os testes passaram." \
		|| echo "$$failed teste(s) falharam."; \
	exit $$failed

# Regra genérica: compila tests/foo.cpp linkando contra a biblioteca com Google Test
build/tests/%: tests/%.cpp $(LIB)
	@mkdir -p build/tests
	$(CXX) $(TEST_CPPFLAGS) $< -L. -lls $(TEST_LDFLAGS) -o $@

# Regra especial: Octave comparison tests.
# octave_ref.cpp (única TU que inclui os cabeçalhos do Octave) é compilada junto.
TEST_OCTAVE_REF := tests/octave_ref.cpp
$(TEST_OCTAVE_BIN): $(TEST_OCTAVE_SRC) $(TEST_OCTAVE_REF) tests/octave_ref.hpp $(LIB)
	@mkdir -p build/tests
	@if [ -z "$(OCTAVE_OCTLIBDIR)" ]; then \
		echo "ERROR: Octave nao encontrado. Instale o pacote octave-dev (ou liboctave-dev)."; \
		exit 1; \
	fi
	$(CXX) $(TEST_CPPFLAGS) $(TEST_OCTAVE_SRC) $(TEST_OCTAVE_REF) -L. -lls $(TEST_LDFLAGS) -o $@

-include $(DEPS)

# --- Benchmark ---
BENCHMARK_BIN := build/benchmark

benchmark: $(BENCHMARK_BIN)
	./$(BENCHMARK_BIN)

$(BENCHMARK_BIN): tests/benchmark.cpp $(LIB)
	@mkdir -p build
	$(CXX) $(CFLAGS) $< -L. -lls $(LDFLAGS) -o $@

<<<<<<< HEAD
# --- Benchmark Runner (TCC: gera CSVs em data/ para os gráficos em Python) ---
BENCHMARK_RUNNER_BIN := build/benchmark_runner

# Compila e executa o runner a partir da raiz (necessário para o caminho data/).
benchmark-runner: $(BENCHMARK_RUNNER_BIN)
	@mkdir -p data
	./$(BENCHMARK_RUNNER_BIN)

$(BENCHMARK_RUNNER_BIN): tests/benchmark_runner.cpp $(LIB)
	@mkdir -p build
	$(CXX) $(CFLAGS) $< -L. -lls $(LDFLAGS) -o $@
=======
COMPARE_CSV    := plot/comparar.csv
COMPARE_PNG    := plot/comparar.png
ART_CSV        := plot/art_data.csv
ARTIGO_IMG_DIR := artigo2/imagens
ARTIGO_TEMPO   := $(ARTIGO_IMG_DIR)/grafico_tempo_execucao.png
ARTIGO_RESIDUO := $(ARTIGO_IMG_DIR)/grafico_norma_residuo.png

# Compila e roda as variantes blas e lapack, mescla os CSVs e gera o gráfico
comparar:
	$(MAKE) clean && $(MAKE) $(BENCHMARK_BIN) BLAS=1 LAPACK=0 VARIANTE=blas
	./$(BENCHMARK_BIN) > $(COMPARE_CSV)
	$(MAKE) clean && $(MAKE) $(BENCHMARK_BIN) BLAS=1 LAPACK=1 VARIANTE=lapack
	./$(BENCHMARK_BIN) | tail -n +2 >> $(COMPARE_CSV)
	$(PYTHON) plot/plot.py $(COMPARE_CSV) $(COMPARE_PNG) --split

# Filtra o CSV completo para as séries usadas no artigo
$(ART_CSV): $(COMPARE_CSV)
	$(PYTHON) -c "\
import pandas as pd; \
df = pd.read_csv('$<'); \
mask = ((df['metodo']=='lu') & (df['variante']=='blas')) \
     | ((df['metodo']=='lapack_lu') & (df['variante']=='lapack')); \
df[mask].to_csv('$@', index=False)"

# Gera as figuras do artigo diretamente em artigo2/imagens/
$(ARTIGO_TEMPO): $(ART_CSV) plot/plot.py
	$(PYTHON) plot/plot.py $(ART_CSV) plot/_artigo.png --split
	mv plot/_artigo_tempo.png $(ARTIGO_TEMPO)
	mv plot/_artigo_residuo.png $(ARTIGO_RESIDUO)

figuras: $(ARTIGO_TEMPO)
>>>>>>> main

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

.PHONY: all clean test test-octave benchmark benchmark-runner install uninstall
