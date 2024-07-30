#cd include/WFA2-lib; make clean setup lib_wfa; cd ../../

WFADIR=include/WFA2-lib
HCLUST=include/hclust-cpp

INC = -L$(WFADIR)/lib -I$(WFADIR) -L${HCLUST} -I${HCLUST}
CFLAGS=-lz -lm -fopenmp -lwfacpp
CC_FLAGS=-g -Wall -O3
LDFLAGS=$(CFLAGS)
CXX=g++
CC=gcc
CXX_FLAGS=-std=c++17 -Wall -O2 -g

ODIR = build
SDIR = src

OUT = otter

CSRCFILES := $(wildcard src/*.c)
CPPSRCFILES := $(wildcard src/*.cpp)
CPPSRCFILES_HCLUST := $(wildcard $(HCLUST)/*.cpp)
SRCFILES := $(CSRCFILES) $(CPPSRCFILES) $(CPPSRCFILES_HCLUST)

COBJS= $(patsubst $(SDIR)/%.c,$(ODIR)/%.o,$(CSRCFILES))
CPPOBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(CPPSRCFILES))
CPPOBJS += $(patsubst $(HCLUST)/%.cpp,$(ODIR)/%.o,$(CPPSRCFILES_HCLUST))
OBJS= $(COBJS) $(CPPOBJS)

DEPS = $(OBJS:%.o=%.d)

.PHONY: all clean

all: $(ODIR)/$(OUT)

clean:
	$(RM) $(ODIR)/$(OUT) $(OBJS) $(DEPS)

$(ODIR)/$(OUT): $(OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(INC) $(CFLAGS)

-include $(DEPS)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CC_FLAGS) -MMD -MP -c -o $@ $<

$(ODIR)/%.o: $(HCLUST)/%.cpp
	$(CXX) $(CXX_FLAGS) -MMD -MP -c -o $@ $< $(INC) $(LDFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CXX_FLAGS) -MMD -MP -c -o $@ $< $(INC) $(LDFLAGS)