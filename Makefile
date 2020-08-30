ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

CXX   = mpicxx
CFLAG = -std=c++11 -lmpi
#CXX   = g++
#CFLAG = -std=c++11
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
TOBJ  = $(ODIR)/main.o $(ODIR)/mat2.o
POBJ  = $(ODIR)/plot_v.o $(ODIR)/mat2.o

tully1.x : $(TOBJ)
	$(CXX) -o $@ $^ $(CFLAG)
pot.x : $(POBJ)
	$(CXX) -o $@ $^ $(CFLAG)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean
clean:
	rm -rf *.x $(ODIR)
