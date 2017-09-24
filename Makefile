GDAL_CFLAGS ?= $(shell gdal-config --cflags)
GDAL_LDFLAGS ?= $(shell gdal-config --libs) $(shell gdal-config --dep-libs)
CC = gcc
CFLAGS ?= -Wall -march=native -mtune=native -Ofast -g
LDFLAGS += -lproj $(GDAL_LDFLAGS)


all: moo

moo: main.o load.o pngwrite.o
	$(CC) $(LDFLAGS) main.o load.o pngwrite.o -o $@

load.o: load.c
	$(CC) $(CFLAGS) $(GDAL_CFLAGS) $< -c -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f *.o

cleaner: clean
	rm -f moo

cleanest: cleaner
