CXXFLAGS    = -g -std=c++11 -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -fdiagnostics-color=always

CPPFLAGS = `root-config --cflags --ldflags`
LDFLAGS = `root-config --libs`

all: deapfit deaplib

deapfit: ./DEAP_fitting.cpp DEAPFitFunction.cpp DEAPFitFunction.h sarge.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ./DEAP_fitting.cpp DEAPFitFunction.cpp sarge.cpp -o deapfit

makeplots: ./make_plots.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ./make_plots.cpp -o makeplots

makegains: ./make_hv_vs_gains.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ./make_hv_vs_gains.cpp -o makegains

deaplib: ./DEAP_fitting.cpp DEAPFitFunction.cpp DEAPFitFunction.h DEAPFitFunction_RootDict.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -shared -fPIC $(LDFLAGS) $^ -o libdeapfit.so

DEAPFitFunction_RootDict.cpp: DEAPFitFunction_Linkdef.h
	@echo "making $@"
	rootcint -f $@ -c -p $(CPPFLAGS) DEAPFitFunction.h $^
