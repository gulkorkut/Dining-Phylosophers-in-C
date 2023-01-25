CC = gcc
CFLAGS = -pthread

all: phsp

phsp:phsp.c 
	$(CC) $(CFLAGS) -o phsp phsp.c -lm -lpthread

clean:
	rm -f phsp

run:
	./phsp 5 500 1000 50 100 exponential 100