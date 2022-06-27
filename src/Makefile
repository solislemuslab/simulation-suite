main: main.cpp SimSuite.cpp SimSuite.hpp core/Network.cpp core/Network.hpp core/Node.cpp core/Node.hpp core/MSEvents.hpp
	g++ -Wall -std=c++11 -g SimSuite.cpp SimSuite.hpp core/MSEvents.hpp main.cpp core/Network.cpp core/Network.hpp core/Node.cpp core/Node.hpp -o main

core-main: core/core-main.cpp SimSuite.cpp SimSuite.hpp core/Network.cpp core/Network.hpp core/Node.cpp core/Node.hpp core/MSEvents.hpp
	g++ -Wall -std=c++11 -g SimSuite.cpp SimSuite.hpp core/MSEvents.hpp core/core-main.cpp core/Network.cpp core/Network.hpp core/Node.cpp core/Node.hpp -o core-main