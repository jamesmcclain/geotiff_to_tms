CC = gcc
CXX = g++
BOOST_CXXFLAGS ?= -I${HOME}/local/src/boost_1_65_1
GDAL_CFLAGS ?= $(shell gdal-config --cflags)
GDAL_LDFLAGS ?= $(shell gdal-config --libs) $(shell gdal-config --dep-libs)
CFLAGS ?= -Wall -Wextra -O0 -ggdb3
CXXFLAGS ?= -std=c++17 $(CFLAGS)
LDFLAGS += -lproj $(GDAL_LDFLAGS)


all: landsat-index landsat-server

landsat-index: landsat-index.o projection.o
	$(CXX) $^ $(LDFLAGS) -fopenmp -o $@

landsat-server: server.o landsat-server.o pngwrite.o fullio.o projection.o
	$(CXX) $^ $(LDFLAGS) -fopenmp -o $@

landsat-server.o: landsat-server.cpp load.h constants.h landsat_scene_handles.hpp greater_landsat_scene.h lesser_landsat_scene.h rtree.hpp textures.hpp
	$(CXX) -fopenmp  $(CXXFLAGS) -Wno-reorder -Wno-unused-parameter $(GDAL_CFLAGS) $(BOOST_CXXFLAGS) $< -c -o $@

landsat-index.o: landsat-index.cpp lesser_landsat_scene.h rtree.hpp
	$(CXX) -fopenmp $(CXXFLAGS) $(GDAL_CFLAGS) $(BOOST_CXXFLAGS) $< -c -o $@

projection.o: projection.c projection.h
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
