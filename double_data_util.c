#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include <glib.h>
#include "config.h"
#include "gts.h"
#include "double_data_util.h"


void read_gdouble_array(
     gdouble   **p_phi,
     gint       *n,
     gchar      *fname,
     gint     num_cells_to_peel)
{
   FILE *fp;
   gint  nxyz, nxy, i, j, k;
   gdouble *phi, *phi1;
   
   gchar     *fbase,command[256];
   gboolean  zip, file_exists;
   gint    n_new[3], nxy_new, nxyz_new, ind, ind1;
   
   file_exists = g_file_test(fname,(GFileTest)G_FILE_TEST_EXISTS);

   if( file_exists )
   {
       zip = g_str_has_suffix(fname, ".gz");       //test if file is zipped

       if( zip )
       { 
           sprintf(command,"gunzip -f %s",fname);
	   system(command);
	   fbase = g_strndup(fname, strlen(fname)-3);
       }
       else fbase = g_strndup(fname, strlen(fname));
   }
   else
   {
       g_error("\nFile %s doesn't exist.",fname);
       return;
   }
   
   fp = fopen(fbase,"r");
   fread(n,sizeof(gint),3,fp); //read dimensions n[0],n[1],n[2]
   
   nxy =  n[0]*n[1];
   nxyz = nxy*n[2];
   phi = (gdouble *)g_malloc(nxyz*sizeof(gdouble));
   fread(phi,sizeof(gdouble),nxyz,fp);
   fclose(fp);
   
   if(num_cells_to_peel)
   {            
      n_new[0] = n[0] - 2*num_cells_to_peel;
      n_new[1] = n[1] - 2*num_cells_to_peel;
      n_new[2] = n[2] - 2*num_cells_to_peel;
      nxy_new = n_new[0] * n_new[1];
      nxyz_new = nxy_new * n_new[2];
      
      phi1 = (gdouble *)g_malloc(nxyz_new*sizeof(gdouble));
      ind1 = 0;
      for(k = num_cells_to_peel; k < n[2] - num_cells_to_peel; k++)
      {
        for(j = num_cells_to_peel; j < n[1] - num_cells_to_peel; j++)
        {
           for(i = num_cells_to_peel; i < n[0] - num_cells_to_peel;i++,ind1++)
	   {
	      ind = i + j*n[0] + k*nxy;
	      phi1[ind1] = phi[ind];
	   }
        }
          
      }

      for(i=0; i < 3; i++) n[i] = n_new[i]; nxyz = nxyz_new;
      g_free(phi);
      phi = phi1;
      printf("\nPeeling done.");
      printf("\nnxyz %d, %d %d %d\n",nxyz,n[0],n[1],n[2]);
  }    
  
  *p_phi = phi;
}


void write_gdouble_array(gdouble *phi,gint *n,gchar *fname)
{
   FILE *fp;
   gint  nxyz;
   
   fp = fopen(fname,"w");
   fwrite(n,sizeof(gint),3,fp); //read dimensions n[0],n[1],n[2]
   
   nxyz = n[0]*n[1]*n[2];
   fwrite(phi,sizeof(gdouble),nxyz,fp);
   fclose(fp); 
}


void read_gfloat_array_into_gdouble(
     gdouble **p_phi,
     gint     *n,
     gchar    *fname,
     gint     num_cells_to_peel)
{
   FILE    *fp;
   gint     nxyz, nxy, i, j, k;
   gfloat  *fphi;
   gdouble *phi, *phi1;
   gint    n_new[3], nxy_new, nxyz_new, ind, ind1;
   
   gchar     *fbase,command[256];
   gboolean  zip, file_exists;
   
   guint   pad = 0;
   gdouble pad_value, temp;
   
   file_exists = g_file_test(fname,(GFileTest)G_FILE_TEST_EXISTS);

   if( file_exists )
   {
       zip = g_str_has_suffix(fname, ".gz");       //test if file is zipped

       if( zip )
       { 
           sprintf(command,"gunzip -f %s",fname);
	   system(command);
	   fbase = g_strndup(fname, strlen(fname)-3);
       }
       else fbase = g_strndup(fname, strlen(fname));
   }
   else
   {
       g_error("\nFile %s doesn't exist.",fname);
       return;
   }
   
   fp = fopen(fbase,"r");
   fread(n,sizeof(gint),3,fp); //read dimensions n[0],n[1],n[2]
   
   nxy =  n[0]*n[1];
   nxyz = nxy * n[2];
   fphi = (gfloat *)g_malloc(nxyz*sizeof(gfloat));
   fread(fphi,sizeof(gfloat),nxyz,fp);
   fclose(fp);
   if( zip )
   { 
     sprintf(command,"gzip -f %s",fbase);
     system(command);
   }
   
   
   phi = (gdouble *)g_malloc(nxyz*sizeof(gdouble));
   for(i = 0; i < nxyz; i++) phi[i] = (gdouble)fphi[i];
   g_free(fphi);
   //printf("\nnxyz %d, %d %d %d\n",nxyz,n[0],n[1],n[2]);
   //fflush(stdout); fflush(stderr);
      
   if( num_cells_to_peel)
   {            
      n_new[0] = n[0] - 2*num_cells_to_peel;
      n_new[1] = n[1] - 2*num_cells_to_peel;
      n_new[2] = n[2] - 2*num_cells_to_peel;
      nxy_new = n_new[0] * n_new[1];
      nxyz_new = nxy_new * n_new[2];
      
      phi1 = (gdouble *)g_malloc(nxyz_new*sizeof(gdouble));
      ind1 = 0;
      for(k = num_cells_to_peel; k < n[2] - num_cells_to_peel; k++)
      {
        for(j = num_cells_to_peel; j < n[1] - num_cells_to_peel; j++)
        {
           for(i = num_cells_to_peel; i < n[0] - num_cells_to_peel;i++,ind1++)
	   {
	      ind = i + j*n[0] + k*nxy;
	      phi1[ind1] = phi[ind];
	   }
        }
          
      }

      for(i=0; i < 3; i++) n[i] = n_new[i]; nxyz = nxyz_new;
      g_free(phi);
      phi = phi1;
      printf("\nPeeling done.");
      printf("\nnxyz %d, %d %d %d\n",nxyz,n[0],n[1],n[2]);
  }
  
   
  if (pad)
  { 
  
    pad_value = 0.02;
  
    /*for(i = 0; i < nxyz; i++)
    {
     temp = fabs(phi[i]);
    
     if( temp > 0 && temp < pad_value) pad_value = temp;
    }*/
    ind = 0;
    for(k = 0; k < n[2]; k++)
      {
        for(j = 0; j < n[1]; j++)
        {
           for(i = 0; i < n[0];i++,ind++)
	   {
	      
	      if( (i == 0)    || (j == 0)    || (k == 0)     ||
	          (i == n[0]) || (j == n[1]) || (k == n[2]) )
	      phi[ind] = pad_value;
	   }
        }
          
      }
  
    printf("\nPadding with value %g done.", pad_value);
  }
  
  
  *p_phi = phi;
}
			
