CC=cc
CXX=c++
CFLAGS=-std=c11 -Wall -O3
CXXFLAGS=-std=c++11 -Wall -O3
CXX_LDFLAGS=-lgsl -lgslcblas
CC_LDFLAGS=-lm
all : gen_fits gen_avgs

gen_fits : fit.o
	$(CXX) $(CXXFLAGS) -o gen_fits fit.o $(CXX_LDFLAGS)
fit.o : fit.cpp minimizer.hpp
	$(CXX) $(CXXFLAGS) -c fit.cpp

gen_avgs : parse_csv.o
	$(CC) $(CFLAGS) -o gen_avgs parse_csv.o $(CC_LDFLAGS)
parse_csv.o : parse_csv.c
	$(CC) $(CFLAGS) -c parse_csv.c

clean :
	rm -f gen_fits gen_avgs fit.o parse_csv.o
