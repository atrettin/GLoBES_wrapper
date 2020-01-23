#connecting globes
prefix = /afs/ifh.de/user/t/trettin/scratch/software/GLoBES_v401/installed/
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

# make sure that lboost_python36 can be found and will be found again
# when the module is loaded without first loading cvmfs 
RLIB=/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/
RLIB64=/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib64/
#BOOST_INC = $(I3_PORTS)/../include/boost/
#BOOST_LIB = $(I3_PORTS)/../lib/boost/

wrapper.so: wrapper.o 
	g++ -shared -Wl,--export-dynamic wrapper.o -lboost_python36 $(PYTHON_INCLUDE) $(PYTHON_LIB)  $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -Wl,-rpath,$(RLIB) -L$(RLIB) -L./. -lsnu -o wrapper.so
	
wrapper.o: wrapper.cxx
	g++ $(PYTHON_INCLUDE) $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -Wl,-rpath,$(RLIB) -fPIC -c wrapper.cxx 

install: wrapper.so
	mv wrapper.so GLoBES.so
	 	 			 
clean:
	rm GLoBES.so wrapper.o	 