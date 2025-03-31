all: main

main: main.o energy.o
	g++ -o $@ $^