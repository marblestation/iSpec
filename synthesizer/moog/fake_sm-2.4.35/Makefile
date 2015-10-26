
all: lib/libfakesm.a


src/fakesm.o: src/fakesm.c
	gcc -c -o src/fakesm.o src/fakesm.c

lib/libfakesm.a: src/fakesm.o
	# This creates the static library. r means to insert with replacement, c means to create a new archive, and s means to write an index.
	mkdir -p lib/
	ar rcs lib/libfakesm.a src/fakesm.o

clean:
	rm -f src/fakesm.o
	rm -f lib/libfakesm.a


