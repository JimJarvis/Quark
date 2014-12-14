EIGEN_PATH = E:\\Dropbox\\Programming\\C\\Eigen-3.2.1
EIGEN_PATH = ../QUARK/lib

CXX = g++
CC = $(CXX)

CXXFLAGS = -std=c++11 -I $(EIGEN_PATH) -O3 -fno-rtti -flto

# detect Cygwin or Unix
ifeq ($(OS),Windows_NT)
	CXXFLAGS += -fwhole-program
else
	CXXFLAGS += -fPIC
endif

LDFLAGS = $(CXXFLAGS)
CFLAGS = $(CXXFLAGS)

QUHEAD = qureg.h qugate.h qumat.h quarklang.h
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

libquark.dylib: $(QUOBJ)
	g++ -dynamiclib -o libquark.dylib $(QUOBJ)

libquark.so: $(QUOBJ)
	g++ -shared -o libquark.so $(QUOBJ)

libquark.dll: $(QUOBJ)
	g++ -shared -o libquark.dll $(QUOBJ)

.PHONY: staticlib
staticlib: libquark.a
	@echo "Static library created"

.PHONY: maclib
maclib: libquark.dylib
	@echo "Shared library on Mac (.dylib) created"

.PHONY: linuxlib
linuxlib: libquark.so
	@echo "Shared library on Linux (.so) created"

.PHONY: winlib
winlib: libquark.dll
	@echo "Shared library on Win (.dll) created"

.PHONY: clean
clean:
	@rm -rf *~ *.o *.exe *.a *.so *.dylib *.dll shor grover basics
	@echo "cleaned"

.PHONY: all
all: shor grover basics
	@echo
	@echo "Recompiled"

# copy headers to quarklang
.PHONY: sync
sync: 
	@echo Sync headers with quarklang 
	@cp -f *.h ../QUARK/lib
