CXXFLAGS    = -g -std=c++11 -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -fdiagnostics-color=always

CPPFLAGS = `root-config --cflags --ldflags`
LDFLAGS = `root-config --libs`

all: deapfit

deapfit: ./DEAP_fitting.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ./DEAP_fitting.cpp -o deapfit
