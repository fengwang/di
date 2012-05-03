####### Compiler, tools and options
CXX           = g++
DEFINES       = -Wall  
INCPATH       = -Iinclude  -I/home/feng/include
LINK          = g++
LFLAGS        = 
DEL_FILE      = rm -f
DEL_DIR       = rmdir
MOVE          = mv -f
MAKE_DIR      = mkdir

####### Output directory
OBJECTS_DIR   = ./obj
BIN_DIR       = ./bin
CHECK_DIR     = ./test

all: grammar_check 

clean: 
	rm -rf $(OBJECTS_DIR)/*
	rm -rf $(BIN_DIR)/*

grammar_check: $(CHECK_DIR)/grammar_check.cc
	$(CXX) -c $(CXXFLAGS) $(DEFINES) $(INCPATH) -o $(OBJECTS_DIR)/grammar_check.o $(CHECK_DIR)/grammar_check.cc
	$(LINK) $(LFLAGS) -o $(BIN_DIR)/grammar_check $(OBJECTS_DIR)/grammar_check.o
