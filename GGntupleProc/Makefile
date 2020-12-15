objects = main.o 
Run3.exe : $(objects)
	g++ -o Run3.exe $(objects) `root-config --libs`
main.o : main.C Run3.h Run3.C
	g++ -Wall -Wno-deprecated `root-config --cflags` main.C -c
clean :
	rm Run3.exe $(objects)
