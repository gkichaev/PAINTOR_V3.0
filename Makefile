CC = g++
OPTS = -std=c++11 -O

curr = "$(PWD)"

all: PAINTOR

PAINTOR: main.cpp
	$(CC) $(OPTS) main.cpp Functions_CoreModel.cpp Functions_Optimize.cpp Functions_IO.cpp -lm -I/$(curr)/eigen/Eigen -o PAINTOR
