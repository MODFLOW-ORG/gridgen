
#-----------------------------------------------------------------------------
# Makefile
#-----------------------------------------------------------------------------

TARGET  = gridgen

#-----------------------------------------------------------------------------
# On Linux,    this should be 'Linux'
# On Mac OS X, this should be 'MAC_OS_X'
#-----------------------------------------------------------------------------

PLATFORM := $(shell uname -s)
MACHINE := $(shell uname -m)

ifeq ($(PLATFORM),Linux)
  ifeq ($(MACHINE),x86_64)
    PLATFORM := LINUX64
  else
    PLATFORM := LINUX
  endif
else
  ifeq ($(PLATFORM),Darwin)
    PLATFORM := MAC_OS_X
  else
    PLATFORM := OTHER
  endif
endif


#-----------------------------------------------------------------------------
# Choose a compiler & its options
#--------------------------------------------------------------------------

CXX       = g++
LINK      = g++
LINKFLAGS = -shared

OPTS = -MMD -W -Wno-deprecated -fPIC  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE=1


#--------------------------------------------------------------------
# Put all together
#--------------------------------------------------------------------
SUBDIRS =  basic intersection grid model packages gui tools tools/shapefile tools/shapefile/shapelib tools/VPass
INCLUDE = $(addprefix -I,$(SUBDIRS))
SRCS=$(wildcard $(addsuffix /*.cpp,$(SUBDIRS))) #cusg_wrap.cpp
OBJS=$(filter-out cusg_wrap.o, ${SRCS:.cpp=.o})

#LIB = $(X_LIB)
GRID_HEADERS = $(wildcard $(addsuffix /*.h, grid))

CFLAGS   = $(OPTS) $(INCLUDE) 
CXXFLAGS = $(CFLAGS) 

.PHONY: release

all: release

release: CXXFLAGS += -O3
release: gridgen

debug: CXXFLAGS += -g 
debug: gridgen


#--------------------------------------------------------------------

#--------------------------------------------------------------------
$(TARGET): $(OBJS) #cusg_wrap.o
	${LINK} ${LINKFLAGS} -o $@ $(OBJS)


gridgen: $(TARGET) gridgenusg.o
	${CXX} -o $@ ${CXXFLAGS} $(OBJS) $(LIB) gridgenusg.o


clean:
	-rm -f $(OBJS) $(TARGET) Dependencies gridgenusg.o gridgenusg


#--------------------------------------------------------------------
.SUFFIXES: .cpp

.cpp.o:
	${CXX} ${CXXFLAGS} -fPIC -c $< -o $@
	cat $*.d >> Dependencies
	rm -f $*.d 

Dependencies:
	touch Dependencies



