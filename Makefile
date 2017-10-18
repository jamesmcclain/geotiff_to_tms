GDAL_CFLAGS ?= $(shell gdal-config --cflags)
GDAL_LDFLAGS ?= $(shell gdal-config --libs) $(shell gdal-config --dep-libs)
CC = gcc
CFLAGS ?= -Wall -Wextra -O0 -ggdb3
LDFLAGS += -lproj $(GDAL_LDFLAGS)


all: landsat-server

landsat-server: server.o landsat.o pngwrite.o fullio.o
	$(CC) $(LDFLAGS) -fopenmp server.o landsat.o pngwrite.o fullio.o -o $@

landsat.o: landsat.c landsat.h load.h constants.h
	$(CC) -fopenmp $(CFLAGS) $(GDAL_CFLAGS) $< -c -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f *.o

cleaner: clean
	rm -f landsat-server

cleanest: cleaner
