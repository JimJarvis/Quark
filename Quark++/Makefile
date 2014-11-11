CXX = g++
CC = $(CXX)

CXXFLAGS = -std=c++11 -I E:\\Dropbox\\Programming\\C\\Eigen-3.2.1 -O3 -fno-rtti -flto -fwhole-program
LDFLAGS = $(CXXFLAGS)
CFLAGS = $(CXXFLAGS)

QUHEAD = qureg.h qugate.h qumat.h

main: main.o algor.o qugate.o qureg.o qumat.o

main.o: algor.h $(QUHEAD)

algor.o: algor.h $(QUHEAD)

qugate.o: $(QUHEAD)

qumate.o: $(QUHEAD)

qureg.o: utils.h

.PHONY: clean
clean:
	@rm -rf *~ *.o Quark
	@echo "cleaned"

.PHONY: all
all: clean main
	@echo
	@echo "Recompiled"
