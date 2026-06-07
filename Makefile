CC      = gcc
CFLAGS  = -Wall -Wextra -std=c11 -g
LDFLAGS = -lm

all: gtfparse

gtfparse: main.o gtf.o
	$(CC) $(CFLAGS) -o gtfparse main.o gtf.o $(LDFLAGS)

main.o: main.c gtf.h
	$(CC) $(CFLAGS) -c main.c

gtf.o: gtf.c gtf.h
	$(CC) $(CFLAGS) -c gtf.c

clean:
	rm -f *.o gtfparse
