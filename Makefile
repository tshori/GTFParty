CC      = gcc
CFLAGS  = -Wall -Wextra -std=c11 -g
LDFLAGS = -lm

all: gtfparse

gtfparse: main.o gtf.o stats.o compare.o
	$(CC) $(CFLAGS) -o gtfparse main.o gtf.o stats.o compare.o $(LDFLAGS)

main.o: main.c gtf.h stats.h compare.h
	$(CC) $(CFLAGS) -c main.c

gtf.o: gtf.c gtf.h
	$(CC) $(CFLAGS) -c gtf.c

stats.o: stats.c gtf.h stats.h
	$(CC) $(CFLAGS) -c stats.c

compare.o: compare.c gtf.h compare.h
	$(CC) $(CFLAGS) -c compare.c

clean:
	rm -f *.o gtfparse
