GDAL_CFLAGS ?= $(shell gdal-config --cflags)
GDAL_LDFLAGS ?= $(shell gdal-config --libs) $(shell gdal-config --dep-libs)
CC = gcc
CXX = g++
CFLAGS ?= -Wall -Wextra -O0 -ggdb3
CXXFLAGS ?= -std=c++14 $(CFLAGS)
LDFLAGS += -lproj $(GDAL_LDFLAGS)


all: landsat-index landsat-server

landsat-index: landsat-index.o projection.o
	$(CXX) $^ $(LDFLAGS) -fopenmp -o $@

landsat-server: server.o landsat-server.o pngwrite.o fullio.o
	$(CC) $^ $(LDFLAGS) -fopenmp -o $@

landsat-server.o: landsat-server.c landsat.h load.h constants.h
	$(CC) -fopenmp $(CFLAGS) $(GDAL_CFLAGS) $< -c -o $@

landsat-index.o: landsat-index.cpp landsat_scene.h
	$(CXX) -fopenmp $(CXXFLAGS) $(GDAL_CFLAGS) $< -c -o $@

projection.o: projection.c projection.h landsat_scene.h
	$(CC) $(CFLAGS) $(GDAL_CFLAGS) $< -c -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f *.o

cleaner: clean
	rm -f landat-index landsat-server

cleanest: cleaner
