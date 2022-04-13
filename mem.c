#include <stdio.h>
#include <glib.h>
#include "mem.h"

/* 
   2022/04
   
   All of the functions in this file have become deprecated.
   Start removing the use of gMemVTable at this time by creating empty     
   functions.
   Search for "DEPRECATE" to see what got turned off. 
   Right now the functions below are somewhat useless wrappers for glib      
   functions.
   

   This file defines MALLLOC, FREE, CALLOC & REALLOC to be used for consistent
   memory profiling.
   Build on existent Glib memory profiling functions and all it does is making
   them visible to all files in the project.
   
   
   
*/

/* DEPRECATE
static GMemVTable	*mem_table;
*/

/* this has to be performed before any other Glib function is used */
void MEM_START(void)
{
// DEPRECATE     mem_table = glib_mem_profiler_table;
// DEPRECATE     g_mem_set_vtable (mem_table);
}

void MEM_STATS(void)
{
// DEPRECATE     g_mem_profile();
}

gpointer MALLOC(gulong	size)
{
        gchar	*new;

// DEPRECATE	new = (gchar *)mem_table->malloc(size);
	new = (gchar *)g_malloc(size);
	
	return new;
}

gpointer CALLOC(gulong size, gulong size_item)
{
        gchar	*new;

//DEPRECATE	new = (gchar *)mem_table->calloc(size,size_item);
	new = (gchar *)g_malloc0_n(size,size_item);
	
	return new;
}

gpointer REALLOC(
	gpointer ptr,
	gulong	size)
{
        gchar	*new;

// DEPRECATE	new = (gchar *)mem_table->realloc(ptr,size);
	new = (gchar *)g_realloc(ptr,size);
	
	return new;
}

void FREE(gpointer ptr)
{
// DEPRECATE        mem_table->free(ptr);
        g_free(ptr);
}
