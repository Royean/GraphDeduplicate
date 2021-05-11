CXX = g++
CXXFLAGS = -Wall

cpres:
	$(CXX) $(CXXFLAGS) -std=c++11 -O3 -o cpres3 GraphDeduplicator.cpp Graph.cpp main.cpp
	
clean:
	rm cpres3
