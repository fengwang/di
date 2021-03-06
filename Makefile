####### Compiler, tools and options
#CXX           = g++
#DEFINES       = -w -std=c++11 -pg
#CXXFLAGS      = 
INCPATH       = -Iinclude  -I/home/feng/include -I/Users/feng/include
LINK          = $(CXX)
#LFLAGS        = -pg
DEL_FILE      = rm -f
DEL_DIR       = rmdir
MOVE          = mv -f
MAKE_DIR      = mkdir

####### Output directory
OBJECTS_DIR   = ./obj
BIN_DIR       = ./bin
CHECK_DIR     = ./test

all: grammar_check ini_test

clean: 
	rm -rf $(OBJECTS_DIR)/*
	rm -rf $(BIN_DIR)/*

grammar_check: $(CHECK_DIR)/grammar_check.cc
	$(CXX) -c $(CXXFLAGS) $(DEFINES) $(INCPATH) -o $(OBJECTS_DIR)/grammar_check.o $(CHECK_DIR)/grammar_check.cc
	$(LINK) $(LFLAGS) -o $(BIN_DIR)/grammar_check $(OBJECTS_DIR)/grammar_check.o

ini_test: $(CHECK_DIR)/test.cc
	$(CXX) -c $(CXXFLAGS) $(DEFINES) $(INCPATH) -o $(OBJECTS_DIR)/ini_test.o $(CHECK_DIR)/test.cc
	$(LINK) $(LFLAGS) -o $(BIN_DIR)/ini_test $(OBJECTS_DIR)/ini_test.o
