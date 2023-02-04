#This code compiles the object files in DEBUGGING mode and links them for execution
g++ -std=c++11 -pthread -lboost_system -lboost_filesystem-mt -g -c *.cpp
g++ -std=c++11 -pthread -lboost_system -lboost_filesystem-mt -g *.o
