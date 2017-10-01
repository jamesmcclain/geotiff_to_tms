GDAL_CFLAGS ?= $(shell gdal-config --cflags)
GDAL_LDFLAGS ?= $(shell gdal-config --libs) $(shell gdal-config --dep-libs)
CC = gcc
CFLAGS ?= -Wall -Werror -O0 -ggdb3
LDFLAGS += -lproj $(GDAL_LDFLAGS)


all: landsat-server

landsat-server: server.o landsat.o pngwrite.o
	$(CC) $(LDFLAGS) server.o landsat.o pngwrite.o -o $@

landsat.o: landsat.c load.h constants.h
	$(CC) $(CFLAGS) $(GDAL_CFLAGS) $< -c -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f *.o

cleaner: clean
	rm -f landsat-server

cleanest: cleaner
