#connecting globes
prefix = /afs/ifh.de/user/t/trettin/scratch/software/GLoBES_v401/installed/
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config
 
local_CFLAGS = -g -O4
 
INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
  
libsnu.o: snu.c
	gcc -std=c99 $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -fPIC -c snu.c -o libsnu.o
  	 
libsnu.so: libsnu.o
	gcc -shared -Wl,--export-dynamic libsnu.o $(INCFLAGS) $(LDFLAGS) $(local_LDFLAGS) -o libsnu.so
   	 	 
install: libsnu.so  libsnu.o snu.c
	cp libsnu.o ${prefix}/lib
	cp libsnu.so ${prefix}/lib
	cp snu.h ${prefix}/include
	 	 			 
clean:
	rm libsnu.so libsnu.o   	 	 			 	 
