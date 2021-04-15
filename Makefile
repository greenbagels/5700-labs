CC=cc
CXX=c++
CFLAGS=-std=c11 -Wall -O3
CXXFLAGS=-std=c++14 -Wall -O3
CXX_LDFLAGS=-lgsl
CC_LDFLAGS=-lm
all : gen_fits gen_avgs

gen_fits : fit.o
	$(CXX) $(CXXFLAGS) $(CXX_LDFLAGS) -o gen_fits fit.o
fit.o : fit.cpp minimizer.hpp
	$(CXX) $(CXXFLAGS) -c fit.cpp

gen_avgs : parse_csv.o
	$(CC) $(CFLAGS) $(CC_LDFLAGS) -o gen_avgs parse_csv.o
parse_csv.o : parse_csv.c
	$(CC) $(CFLAGS) -c parse_csv.c

clean :
	rm -f gen_fits gen_avgs fit.o parse_csv.o
