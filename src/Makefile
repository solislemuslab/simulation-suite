CPPFLAGS := -Wall -std=c++11 -g
CPP := g++

ifeq ($(OS),Windows_NT)
	RM := del /f
endif


build: core-main main

all: core-main main test

test: CoreTests.o SimSuite.o Network.o Node.o
	$(CPP) $(CPPFLAGS) CoreTests.o SimSuite.o Network.o Node.o -o test


Node.o: core/Node.cpp core/Node.hpp
	$(CPP) $(CPPFLAGS) -c core/Node.cpp core/Node.hpp

Network.o: core/Network.cpp core/Network.hpp core/MSEvents.hpp
	$(CPP) $(CPPFLAGS) -c core/Network.cpp core/Network.hpp core/MSEvents.hpp

SimSuite.o: SimSuite.cpp SimSuite.hpp
	$(CPP) $(CPPFLAGS) -c SimSuite.cpp SimSuite.hpp

CoreTests.o: tests/CoreTests.cpp tests/CoreTests.hpp
	$(CPP) $(CPPFLAGS) -c tests/CoreTests.cpp tests/CoreTests.hpp

core-main: core/core-main.cpp Network.o Node.o
	$(CPP) $(CPPFLAGS) core/core-main.cpp Network.o Node.o -o core-main

main: SimSuite.o Network.o Node.o main.cpp
	$(CPP) $(CPPFLAGS) SimSuite.o Network.o Node.o main.cpp -o main

clean:
	$(RM) *.o
	$(RM) *.exe
	$(RM) ./core/*.gch
	$(RM) ./tests/*.gch
	$(RM) ./*.gch