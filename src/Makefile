all: main_exe clean

CFLAGS = -Wall -O0 -lm -Wextra -g

main_exe: constants.o utils.o physics.o main.o
			gcc -o main_exe constants.o utils.o physics.o main.o $(CFLAGS)

constants.o: constants.c constants.h
			gcc -o constants.o -c constants.c $(CFLAGS)

utils.o: utils.c utils.h
			gcc -o utils.o -c utils.c $(CFLAGS)

physics.o: physics.c physics.h
			gcc -o physics.o -c physics.c $(CFLAGS)

main.o: main.c physics.h utils.h constants.h
			gcc -o main.o -c main.c $(CFLAGS)

clean:
			rm -f *.o core

rm_exe:
			rm -f main_exe core

clean_data:
			rm -f ../data/data/*.txt core
			rm -f ../data/potential/*.txt core
			rm -f ../data/density/*.txt  core
