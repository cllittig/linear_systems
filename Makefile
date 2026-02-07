
CXX := g++
CFLAGS := -std=c++17 -Wall -Wextra -Iinclude
LDFLAGS := 

SOURCE := $(shell find src -name '*.cpp')
OBJS := $(SOURCE:src/%.cpp=build/%.o)


#linka o código-objeto
main: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS}

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CFLAGS} -c $< -o $@

clean:
	rm -rf main build
