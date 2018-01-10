CC = g++
OPTS = -std=c++11 -O

curr = "$(PWD)"

all: PAINTOR
PAINTOR: main.cpp
	$(CC) $(OPTS) main.cpp Functions_model.cpp Functions_optimize.cpp Functions_IO.cpp -I/$(curr)/eigen/Eigen -L/${curr}/lib -I/${curr}/include -lm -lnlopt -std=c++11 -o PAINTOR
