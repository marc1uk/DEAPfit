#!/bin/bash
rm -f ./DEAP_fitting.tar.gz
tar -zcf DEAP_fitting.tar.gz DEAPFitFunction.cpp DEAPFitFunction.h DEAPFitFunction_Linkdef.h DEAP_fitting.cpp Makefile sarge.cpp sarge.h
rm /pnfs/annie/persistent/users/moflaher/deapfit/DEAPfit/DEAP_fitting.tar.gz
cp DEAP_fitting.tar.gz /pnfs/annie/persistent/users/moflaher/deapfit/DEAPfit/DEAP_fitting.tar.gz

