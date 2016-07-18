
CC = g++

#CFLAGS = -Wall -c -m64 -std=c++11 -O3 
CFLAGS = -Wall -m64 -std=c++11 -O3 
LFLAGS = -Wall -m64 -std=c++11 

#OBJS = Main.o

#STRUCTUREify: $(OBJS)
#	$(CC) $(LFLAGS) $(OBJS) -o STRUCTUREify

STRUCTUREify: STRUCTUREify.cpp
	$(CC) $(CFLAGS) STRUCTUREify.cpp -o STRUCTUREify

clean:
	rm -rf *.o STRUCTUREify
