#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

/* 2022/04

   Should be removed from new contact angle code installation.
    
   This is memory profiler that used to work, but is now depracated.
   I do not see any options for replacement at this time.

*/

main()
{
     GMemVTable	*mem_table;
     gint *a;
     
     mem_table = glib_mem_profiler_table;
     g_mem_set_vtable (mem_table);
     
     a = mem_table->malloc(100*sizeof(gint));
    
     mem_table->free(a); 
     
     g_mem_profile();
     
}
