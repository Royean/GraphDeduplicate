CXX = g++
CXXFLAGS = -Wall

cpres:
	$(CXX) $(CXXFLAGS) -std=c++11 -O3 -o cpres3  main.cpp GraphDeduplicator.cpp Bi2Tri.cpp Utility.cpp FPGrowth.cpp
	
debug:
	$(CXX) $(CXXFLAGS) -std=c++11 -O3 -o cpres3  main.cpp GraphDeduplicator.cpp Bi2Tri.cpp Utility.cpp -g

clean:
	rm cpres3
