#connecting globes
prefix = /afs/ifh.de/user/t/terliuk/scratch/software/GLoBES_v311/installed//
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)

# location of the Python header files
PYTHON_INCLUDE:=$(shell python-config --includes)
PYTHON_LIB:=$(shell python-config --libs)
#PYTHON_LIB =  $(I3_PORTS)/../lib/python$(PYTHON_VERSION)
# location of the Boost Python include files and library
 
#BOOST_INC = $(I3_PORTS)/../include/boost/
#BOOST_LIB = $(I3_PORTS)/../lib/boost/

TARGET=wrapper
$(TARGET).so: $(TARGET).o 
	g++ -shared -Wl,--export-dynamic  $(TARGET).o -lboost_python $(PYTHON_INCLUDE) $(PYTHON_LIB)  $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -L./. -lsnu -o $(TARGET).so ; \
    mv wrapper.so GLoBES.so

$(TARGET).o: $(TARGET).cxx
	g++ $(PYTHON_INCLUDE) $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -fPIC -c $(TARGET).cxx 

