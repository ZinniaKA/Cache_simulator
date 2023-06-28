CXX=g++
CFLAGS=-I. -std=c++17 -Wall -Werror -Wextra

%.o: %.cpp 
	$(CXX) -c -o $@ $< $(CFLAGS)

all: main.o  
	$(CXX) -o cache_simulate main.o $(CFLAGS)
clean: 
	rm -f main.o cache_simulate