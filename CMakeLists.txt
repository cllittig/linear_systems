cmake_minimum_required(VERSION 3.14)
project(SistemasLineares VERSION 1.0)

# Configurações básicas
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Opções de compilação
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")
endif()

# Inclui diretórios de headers
include_directories(
    include
    include/algebra_linear
    include/metodos
    include/utils
)

# Biblioteca de Álgebra Linear
add_library(algebra_linear
    src/algebra_linear/matriz.cpp
    src/algebra_linear/sistema_linear.cpp
)

# Biblioteca de Métodos Numéricos
add_library(metodos
    src/metodos/gauss_jordan
)

# Biblioteca de Utilitários
add_library(utils
    src/utils/read_data.cpp
)

# Executável principal
add_executable(sistema_linear
    src/main.cpp
)

# Link com as bibliotecas
target_link_libraries(sistema_linear
    algebra_linear
    metodos
    utils
) 

# # Configuração do Google Test (OPCIONAL)
include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  

# # Testes dos Métodos
# # add_executable(test_metodos
# #     tests/unit/metodos/test_gauss.cpp
# #     tests/unit/metodos/test_gauss_seidel.cpp
# # )

# target_link_libraries(test_metodos
#     algebra_linear
#     metodos_numericos
#     utils
#     gtest_main
# )

# # Habilita testing
# enable_testing()
# add_test(NAME test_algebra COMMAND test_algebra)
# add_test(NAME test_metodos COMMAND test_metodos)

# # Instalação (OPCIONAL)
# install(TARGETS sistema_linear DESTINATION bin)
# install(DIRECTORY include/ DESTINATION include)