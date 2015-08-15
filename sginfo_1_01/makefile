CC = cc
CCFLAGS = -O

SGINFOFLAGS =#-DNo_ListRotMxInfo -DUSE_GS_SI#-DDEBUG_Find_si

sginfo: sginfo.o sgclib.o sgio.o sgfind.o sghkl.o sgsi.o
	$(CC) $(CCFLAGS) -o sginfo sginfo.o sgclib.o sgio.o sgfind.o sghkl.o sgsi.o -lm

sginfo.o: sginfo.c sginfo.h
	$(CC) $(CCFLAGS) -c sginfo.c $(SGINFOFLAGS)

sgtest: sgtest.o sgclib.o sgio.o sgfind.o
	$(CC) $(CCFLAGS) -o sgtest sgtest.o sgclib.o sgio.o sgfind.o

sgtest.o: sginfo.h

sgquick: sgquick.o sgclib.o sgio.o
	$(CC) $(CCFLAGS) -o sgquick sgquick.o sgclib.o sgio.o

sgquick: sginfo.h

sgclib.o: sginfo.h
sgio.o:   sginfo.h
sgfind.o: sginfo.h
sghkl.o:  sginfo.h
sgsi.o:   sginfo.h

libsginfo: sgclib.o sgio.o sgfind.o sghkl.o sgsi.o
	rm -f libsginfo.a
	ar r libsginfo.a sgclib.o sgio.o sgfind.o sghkl.o sgsi.o

.c.o:
	$(CC) $(CCFLAGS) -c $*.c

all: sginfo sgtest sgquick libsginfo

clean:
	rm -f sgclib.o sgio.o sgfind.o sghkl.o sgsi.o
	rm -f sginfo sginfo.o
	rm -f sgtest sgtest.o
	rm -f sgquick sgquick.o
	rm -f libsginfo.a
