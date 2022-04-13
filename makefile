CC = gcc -O2
#CFLAGS = -g
#CFLAGS = -m32
#CFLAGS = -g -O2 -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations

# glib-2.0 is used both in GTS and in this code
#GLIB_CFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
#GLIB_LIBS = -lglib-2.0
GLIB_CFLAGS = `pkg-config --cflags glib-2.0`
GLIB_LIBS = `pkg-config --libs glib-2.0`


# original directory of GTS code (for .h files)
GTS_DIR = $(HOME)/local/src/gts-0.7.6/
GTS_DIR_SRC = $(HOME)/local/src/gts-0.7.6/src

# directory where GTS library is located
LIB_DIR = $(HOME)/local/lib
LIB_ADD = $(LIB_DIR)/libgts.a

INCLUDE = -I. -I$(GTS_DIR) -I$(GTS_DIR_SRC) -I$(HOME)/local/include/

# where to install executable
BIN_DIR = $(HOME)/local/bin/

# object list for contact angle
COBJLIST = contact_angle_main.o cont_angle.o double_data_util.o intersect_surf.o \
	   segfl_util.o viz_util.o mem.o

# object list for extrapolation
EOBJLIST = extrapolate_main.o double_data_util.o iso_modified.o \
	intersect_surf.o viz_util.o mem.o

# object list for interpolation
IOBJLIST =interpolate_main.o double_data_util.o intersect_surf.o iso_modified.o \
	 mem.o coarsen_surface.o viz_util.o 

# object list for test_sphere
TOBJLIST = test_sphere.o cont_angle.o intersect_surf.o mem.o \
	   segfl_util.o viz_util.o

#EXECUTABLES 
all:    change_view_all_steps contact_angle extrapolate interpolate interpolate_all_steps_double_data \
	process_all_blobs process_all_steps_double_data test_sphere

change_view_all_steps:
	$(CC) $(CFLAGS) change_view_all_steps.c -o change_view_all_steps
	mv  change_view_all_steps $(BIN_DIR) 
		
contact_angle: $(COBJLIST)
	$(CC) $(CFLAGS)  $(COBJLIST) $(GLIB_LIBS) $(LIB_ADD) -o contact_angle -lm
	mv contact_angle $(BIN_DIR) 

extrapolate: $(EOBJLIST)
	$(CC) $(CFLAGS) $(EOBJLIST) $(GLIB_LIBS) $(LIB_ADD) -o extrapolate -lm
	mv extrapolate $(BIN_DIR) 

interpolate: $(IOBJLIST)
	$(CC) $(CFLAGS) $(IOBJLIST) $(GLIB_LIBS) $(LIB_ADD) -o interpolate -lm
	mv interpolate $(BIN_DIR) 
	
interpolate_all_steps_double_data:
	$(CC) $(CFLAGS) interpolate_all_steps_double_data.c -o interpolate_all_steps_double_data
	mv interpolate_all_steps_double_data $(BIN_DIR)		
		
process_all_blobs:
	$(CC) $(CFLAGS) process_all_blobs.c -o process_all_blobs
	mv process_all_blobs $(BIN_DIR) 

process_all_steps_double_data:
	$(CC) $(CFLAGS) process_all_steps_double_data.c -o process_all_steps_double_data
	mv process_all_steps_double_data $(BIN_DIR) 	

test_sphere: $(TOBJLIST)	
	$(CC) $(CFLAGS) $(TOBJLIST) $(GLIB_LIBS) $(LIB_ADD) -o test_sphere -lm
	mv test_sphere $(BIN_DIR) 

mem_test: mem_test.o
	$(CC) $(CFLAGS) mem_test.o $(GLIB_LIBS) $(LIB_ADD) -o mem_test -lm
		

#AUXILIARY FILES
contact_angle_main.o: contact_angle_main.c segfl_util.h viz_util.h \
	intersect_surf.h cont_angle.h mem.h \
	$(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c contact_angle_main.c

cont_angle.o: cont_angle.c cont_angle.h intersect_surf.h  \
	mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c cont_angle.c

coarsen_surface.o: coarsen_surface.c coarsen_surface.h \
        mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c coarsen_surface.c
	
double_data_util.o: double_data_util.c double_data_util.h mem.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c double_data_util.c

extrapolate_main.o: extrapolate_main.c viz_util.h intersect_surf.h iso_modified.h \
	mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c extrapolate_main.c

interpolate_main.o: interpolate_main.c coarsen_surface.h viz_util.h intersect_surf.h \
	iso_modified.h mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c interpolate_main.c

intersect_surf.o: intersect_surf.c intersect_surf.h  \
	mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c intersect_surf.c

iso_modified.o: iso_modified.c $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c iso_modified.c

	
mem.o: mem.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c mem.c

mem_test.o: mem.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c mem_test.c

segfl_util.o: segfl_util.c segfl_util.h mem.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c segfl_util.c	
	
test_sphere.o: test_sphere.c segfl_util.h intersect_surf.h cont_angle.h \
	viz_util.h mem.h $(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c test_sphere.c
		
viz_util.o: viz_util.c viz_util.h intersect_surf.h  mem.h \
	$(GTS_DIR)/config.h $(GTS_DIR_SRC)/gts.h
	$(CC) $(CFLAGS) $(INCLUDE) $(GLIB_CFLAGS) -c viz_util.c	

clean:
	rm -f *.o
