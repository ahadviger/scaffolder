CC=g++
CFLAGS=-std=c++0x

all: main

main:
	$(CC) $(CFLAGS) main.cpp contig.cpp scaffold.cpp -o main

clean:
	rm main