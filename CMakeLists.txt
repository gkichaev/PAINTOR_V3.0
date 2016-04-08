cmake_minimum_required(VERSION 3.2.2)
project(PAINTOR_3_0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(curr /Users/glebkichaev/Dropbox/PAINTOR_3.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/${curr}/eigen/Eigen -lm -std=c++11")

set(DCMAKE_CXX_COMPILER /usr/bin/g++)

set(DCMAKE_C_COMPILER /usr/bin/gcc)


set(SOURCE_FILES
        main.cpp
        Functions_CoreModel.cpp
        Functions_Optimize.cpp
        Functions_IO.cpp
)


add_executable(PAINTOR_3_0 ${SOURCE_FILES})
