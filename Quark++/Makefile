EIGEN_PATH = E:\\Dropbox\\Programming\\C\\Eigen-3.2.1

CXX = g++
CC = $(CXX)

CXXFLAGS = -std=c++11 -I $(EIGEN_PATH) -O3 -fno-rtti -flto -fwhole-program
LDFLAGS = $(CXXFLAGS)
CFLAGS = $(CXXFLAGS)

QUHEAD = qureg.h qugate.h qumat.h
QUOBJ = qugate.o qureg.o qumat.o

shor: shor.o $(QUOBJ)

grover: grover.o $(QUOBJ)

basics: basics.o $(QUOBJ)

algor.o: algor.h $(QUHEAD)

shor.o: $(QUHEAD)

grover.o: $(QUHEAD)

basics.o: $(QUHEAD)

qugate.o: $(QUHEAD)

qumate.o: $(QUHEAD)

qureg.o: utils.h


libquark.a: $(QUOBJ)
	ar rcs $@ $^

.PHONY: lib
lib: libquark.a
	@echo "Static library created"

.PHONY: clean
clean:
	@rm -rf *~ *.o *.exe shor grover basics
	@echo "cleaned"

.PHONY: all
all: shor grover basics
	@echo
	@echo "Recompiled"
