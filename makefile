all: main

main: main.o energy.o energy_b.o
	g++ -o $@ $^

%_b.c: %.c
	$(HOME)/tapenade/tapenade_3.16/bin/tapenade -b -head 'TotalEnergy()/(x1,v1)' $<
	sed -e '/adStack/ s|^|// |' -i $@
