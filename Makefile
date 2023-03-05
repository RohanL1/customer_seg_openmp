CC = g++-12

default: kmeans

kmeans: kmeans.cpp
	${CC} -O0 -g -Wall -Wextra -Wno-unused-parameter -fopenmp -o $@ kmeans.cpp
clean:
	-rm -vf kmeans
