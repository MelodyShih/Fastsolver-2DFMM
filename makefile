CXX=g++
CXXFLAGS=-std=c++11

HEADERS=boxtree.hpp

%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $(INC) $< -o $@

driver: driver.o boxtree.o 
	$(CXX) driver.o boxtree.o -o driver

