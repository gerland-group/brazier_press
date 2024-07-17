# Cesar L. Pastrana, 2021

CXX := g++
FLAGS := -march=native -std=c++11 -O3 -lm -Wall
OPENMP_FLAG := -fopenmp
DEBUG := -std=c++11 -g -lm -Wall
PROFILING := -pg

BIN_PATH := ./bin
SRC_PATH := ./src
SRCS := $(wildcard ${SRC_PATH}/*.cpp)

#.PHONY: install clean

single: 
	@clear
	@echo -n "Compiling sequential... "
	@${CXX} ${SRCS} ${FLAGS} -o ${BIN_PATH}/brazier
	@echo "Sequential compilation Done!"

parallel: 
	@clear
	@echo -n "Compiling parallel... "
	@${CXX} ${SRCS} ${FLAGS} ${OPENMP_FLAG} -o ${BIN_PATH}/brazier
	@echo "OpenMP-active compilation Done!"

debug:
	@clear
	@echo -n "Compiling in debug mode... "
	@${CXX} ${SRCS} ${DEBUG} -o ${BIN_PATH}/brazier
	@echo "Debug compilation Done!"

profiling:
	@clear
	@echo -n "Compiling in profiling mode (sequential)... "
	@${CXX} ${SRCS} ${FLAGS} ${PROFILING} -o ${BIN_PATH}/brazier
	@echo "Profiling compilation Done!"

profiling_omp:
	@clear
	@echo -n "Compiling in profiling mode (parallel)... "
	@${CXX} ${SRCS} ${FLAGS} ${OPENMP_FLAG} ${PROFILING} -o ${BIN_PATH}/brazier
	@echo "Profiling compilation Done!"


