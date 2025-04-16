all: main

CFLAGS=-O3

main: main.o energy.o energy_b.o energy_b_d.o vtu_write.o
	g++ -o $@ $^

energy_b.c: energy.c
	$(HOME)/tapenade/tapenade_3.16/bin/tapenade -b -head 'TotalEnergy()/(x1,v1)' $<
	sed -e '/adStack/ s|^|// |' -i $@

ode_c.c : energy_b.c ode.c
	cat $^ >$@

ode_c_d.c: ode_c.c
	$(HOME)/tapenade/tapenade_3.16/bin/tapenade -d -head 'ODEstep(res)/(v2)' $<
	sed -e '/adStack/ s|^|// |' -i $@

energy_b_d.c: energy_b.c
	$(HOME)/tapenade/tapenade_3.16/bin/tapenade -d -head 'TotalEnergy_b(x1b,v1b)/(x1,v1)' $<
	sed -e '/adStack/ s|^|// |' -i $@

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<


.PRECIOUS: %.c