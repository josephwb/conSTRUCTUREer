OBJS = main.o utils.o write_results.o info.o

CC = g++

CFLAGS = -Wall -c -std=c++11 -O3
LFLAGS = -Wall -std=c++11

STRUCTUREify: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o STRUCTUREify

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

utils.o: utils.cpp utils.h
	$(CC) $(CFLAGS) utils.cpp

write_results.o: write_results.cpp write_results.h
	$(CC) $(CFLAGS) write_results.cpp

info.o: info.cpp info.h
	$(CC) $(CFLAGS) info.cpp

clean:
	rm -rf *.o STRUCTUREify
