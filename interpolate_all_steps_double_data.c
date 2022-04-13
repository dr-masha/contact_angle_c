#include <stdio.h>
#include <stdlib.h>


/*
   interpolate_all_steps_double_data start_step end_step
   
   Interpolates a sequence of fluid data stored in files data_step%d.gz
   Assumes grain phase is stored in 'mask.c'.

   Requires:    
   'interpolate' - local (this code) executable for interpolation
   ~/local/geomview_templates/template_fracture.list - desired geomview template   
   'rpl' - Unix command for replacing words in text files.
*/   
   
int main(int argc, char **argv)
{
   char command[256], fname[256], fname_gmv[256];

   int i,istart,iend;
   
   if(argc < 3) 
   {
     printf("\nUsage: interpolate_al_steps_double_data start_step end_step\n");
     return 1;
   }
   
   istart = atoi(argv[1]);
   iend = atoi(argv[2]);
   
   system("mkdir geomview");
     
   //sprintf(command,"interpolate data_init mask surf");
  // system(command);
   
  
   for(i=istart; i <= iend; i++)
   {
     sprintf(command,"interpolate data_step%d.gz mask.gz geomview/surf_d",i);
     system(command);
       
     sprintf(command,"mv geomview/surf_d_fluid1.list geomview/surf_d%i.list",i);
     system(command);
     
     sprintf(command,"mv geomview/surf_d_fluid2.list geomview/surf_wet_d%i.list",i);
     system(command);
 
     sprintf(command,"cp ~/local/geomview_templates/template_fracture.list geomview/view_d%d.list",i);
     system(command);
     
     sprintf(command,"rpl 'tmp.list' 'surf_d%i.list' geomview/view_d%d.list",i,i);
     system(command);    

     sprintf(command,"cp ~/local/geomview_templates/template_fracture.list geomview/view_wet_d%d.list",i);
     system(command);

     sprintf(command,"rpl 'tmp.list' 'surf_wet_d%i.list' geomview/view_wet_d%d.list",i,i);
     system(command);    
   }
}
