#include <stdio.h>
#include <stdlib.h>


/*
   change_view_all_steps start_step end_step
   
   Assuming you used 'interpolate_all_steps_double_data' to produce Geomview visualization
   of interpolated data, this program changes Geomview template used for the purpose.
     
   Requires:      
   'rpl' - Unix command for replacing words in text files.
*/   
   
int main(int argc, char **argv)
{
   char command[256];

   int   i, istart, iend;
   
   char  *template_filename;
   
   if(argc < 4) 
   {
     printf("\nUsage: change_view_all_steps start_step end_step new_template\n");
     return 1;
   }
   
   istart = atoi(argv[1]);
   iend = atoi(argv[2]);
   template_filename = argv[3];
          
   for(i=istart; i <= iend; i++)
   {
     sprintf(command,"cp %s geomview/view_d%d.list",template_filename,i);
     system(command);
     
     sprintf(command,"rpl 'tmp.list' 'surf_d%i.list' geomview/view_d%d.list",i,i);
     system(command);    

     sprintf(command,"cp %s geomview/view_wet_d%d.list",template_filename,i);
     system(command);

     sprintf(command,"rpl 'tmp.list' 'surf_wet_d%i.list' geomview/view_wet_d%d.list",i,i);
     system(command);    
   }
}
