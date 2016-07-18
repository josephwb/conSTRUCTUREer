OBJS = main.o utils.o write_results.o

CC = g++

CFLAGS = -Wall -c -m64 -std=c++11 -O3 
LFLAGS = -Wall -m64 -std=c++11 

STRUCTUREify: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o STRUCTUREify

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

utils.o: utils.cpp utils.h
	$(CC) $(CFLAGS) utils.cpp

write_results.o: write_results.cpp write_results.h
	$(CC) $(CFLAGS) write_results.cpp

clean:
	rm -rf *.o STRUCTUREify
