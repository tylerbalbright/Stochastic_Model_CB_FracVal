#This code compiles the object files and links them for execution
g++ -std=c++11 -pthread -lboost_system -lboost_filesystem-mt -c *.cpp
g++ -std=c++11 -pthread -lboost_system -lboost_filesystem-mt *.o
