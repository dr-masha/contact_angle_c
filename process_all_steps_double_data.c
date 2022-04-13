#include <stdio.h>
#include <stdlib.h>

#define CREATE_GMV_FILE(command,fname,fname_gmv,id)\
{\
   sprintf(command,"contact_angle %s %s",fname,color);\
   system(command);\
   \
   sprintf(command,"cp /tmp/surface.list %s",fname_gmv);\
   system(command);\
} 

#define THROAT 0
#define PORE 0 
#define FRACTURE 1


int main(int argc, char **argv)
{
   char command[256], fname[256], fname_gmv[256];

   int i,istart,iend;
   
   if(argc < 3) 
   {
     printf("\nPlease provide starting and ending data_step numbers.");
     return 1;
   }
   
   istart = atoi(argv[1]);
   iend = atoi(argv[2]);
   
   system("mkdir geomview");
     
   sprintf(command,"contact_angle data_init.gz mask.gz surf_init 0");
   system(command);
   
  
   for(i=istart; i <= iend; i++)
   {
     sprintf(command,"contact_angle data_step%d.gz mask.gz geomview/surf_d %d",i,i);
     system(command);
       
     sprintf(command,"cp ~/local/geomview_templates/template_blob1.list geomview/view_d%d.list",i);
     system(command);
    
     
     sprintf(command,"rpl 'tmp.list' surf_d_%i.list geomview/view_d%d.list",i,i);
     system(command);     
   }
}
