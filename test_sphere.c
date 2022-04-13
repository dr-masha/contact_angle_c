#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>

#include "config.h"
#include "gts.h"
#include "intersect_surf.h"
#include "segfl_util.h"
#include "cont_angle.h"
#include "viz_util.h"
#include "mem.h"
  

extern void plot_surface_pack(SurfacePack *,GtsCartesianGrid *,guint32,guint32,
                                                              guint32,gchar *);   
  
// run as 'test_sphere d' where d is diameter, odd integer
int main(int argc, char * argv[])
{

 gint    d;
 gint    i, j, k, ind, nxyz, r, r2, dist1;
 gint    sum1, sum2, sum;
 gfloat  r_eff, r_eff2;
 gfloat  xc, yc, zc, h, hplot;
 
 gfloat  a, b, c, x0, y0, z0, norm;
 
 gdouble mean;
 		
 GtsCartesianGrid g;
 
 guchar *data1, *data2, *data, *data_sphere;
 
 SurfacePack *pack;
 
 guint32  S1 = 1, S2 = 2, SG = 4;
 
 gchar  fname[256],base[256];
 FILE   *fp;
 
 MEM_START();
 
 printf("\n Define a spherical cap of fluid 1 surounded by fluid 2, and in contact with (intersecting) planar grain surface."); 
 printf("\nUsage:");
 printf("\n'test_sphere d'");
 printf("\n'test_sphere d h'");
 printf("\n'test_sphere d h a b c'");
 printf("\n       where d is spherical cap diameter, odd integer");
 printf("\n             h is plotting distance of the interecting plane from sphere center, h < d/2");
 printf("\n             a,b,c define the normal of the intersecting plane");
      
 
 // Sphere diameter
 d = 11;        // default value
 if (argc >= 2) d = atoi(argv[1]);
 printf("\nSphere diameter d = %d",d);
 
 // Distance
 hplot = 0;        // default value
 if (argc >= 3)
 {
    hplot = atof(argv[2]);
 }   
 printf("\nPlotting distance of the intersecting plane from sphere center h = %g",hplot);
 
 // Plane normal
 a = 0; b = 0; c = 1;  //(a,b,c) is default plane normal vector
 if(argc >= 6)
 {
     a = atof(argv[3]);
     b = atof(argv[4]);
     c = atof(argv[5]);
 }
 printf("\nIntersecting plane normal vector (a,b,c)=(%g,%g,%g)",a,b,c);	   
           
 sprintf(base,"sphere%d",d);
 
 /* Initialize grid */
  g.dx = g.dy = g.dz = 1.0; //voxel length assumed 1.0 for now
  g.x = g.y = g.z = -0.5;   //center of the first voxel assumed at (0,0,0)
                            //hence volume boundary starts at (-.5,-.5,-.5)
  g.nx = g.ny = g.nz = d+2; // number of voxels in each dimensions
  
  nxyz = g.nx * g.ny * g.nz; //total number of voxels
  
  data1 = (guchar *)g_malloc(nxyz*UCSZ);
  data2 = (guchar *)g_malloc(nxyz*UCSZ);
  data = (guchar *)g_malloc(nxyz*UCSZ);
  data_sphere = (guchar *)g_malloc(nxyz*UCSZ);
  
  
  /*
  * create digitized data_sphere with diameter 'd'
  * data_sphere == 1 at voxels inside the sphere, otherwise 0
  */
  xc = yc = zc = 0.5*d + 0.5; //sphere center
  printf("\nSphere center (%g,%g,%g)",xc,yc,zc);
  
  if( (xc - floor(xc)) != 0 ) printf("\nMight have issues with diameter, use odd integer.");
   
  r_eff = 0.5*d; r_eff2 = r_eff*r_eff; //effective sphere radius
  r = d/2;  r2 = r*r; //radius that will be used for voxel centers
  ind = 0;
  for(k = 0; k < g.nz; k++) 
    for(j = 0; j < g.ny; j++)
      for(i = 0; i < g.nx; i++, ind++)
      {
          data_sphere[ind] = 0;
	  
	  dist1 = (i - xc)*(i - xc) + (j - yc)*(j - yc) + 
                 (k - zc)*(k - zc);
          if( dist1 < r2 ) data_sphere[ind] = 1;
      }
  
  // File storing distance and mean contact angle computed for this distance.
  sprintf(fname,"sphere%d_theta_vs_dist.csv",d); 
  fp = fopen(fname,"w");
  fprintf(fp,"normalized_distance mean_contact_angle");
  
  // Define intial plane through (xc,yc,zc) and with normal (a,b,c)
  
  norm = sqrt(a*a + b*b + c*c);
  //normalize
  a = a/norm; b = b/norm; c = c/norm;
  
  
  for(h = 0; h < r; h++)
  {
     printf("\n\nDistance h %g",h); fflush(stdout);
     
     //new plane position
     x0 = xc + h*a;
     y0 = yc + h*b;
     z0 = zc + h*c;
      
     ind = 0;
     sum = sum1 = sum2 = 0;
     for(k = 0; k < g.nz; k++) 
      for(j = 0; j < g.ny; j++)
        for(i = 0; i < g.nx; i++, ind++)
        {
          data[ind] = data1[ind] = data2[ind] = 1;
	  
	 
	  if( a*(i - x0) + b*(j - y0) + c*(k - z0) < 0 ) 
	  {
	     data[ind] = 0; //negative side of the plane is considered grain
	     sum++;
	  }
	  else if(data_sphere[ind] == 1)
	  {
	     data1[ind] = 0;  //fluid 1 is inside sphere and above z-level
	     sum1++;
	  }
	  else
	  {
	     data2[ind] = 0;  //fluid 2 is outside sphere and above z-level
	     sum2++;
	  } 
        }

      if( sum && sum1 && sum2)
      {
 
         pack = intersect_surfaces_from_data_base(data1,data2,data,&g,S1,S2,SG,
                                	 FALSE,FALSE);
         mean = process_contact_angles(pack,&g,S1,S2,SG,FALSE,FALSE,base,FALSE,0);
     
         if (h == hplot) plot_surface_pack(pack,&g,S1,S2,SG,base);
      
         destroy_surface_pack(pack);
         fprintf(fp,"\n%g %g",h/(gfloat)r,mean);
       }			 
  }
  fclose(fp);
  printf("\n\n Geomview (*.list) files provide visualization and .csv file theta vs. distance\n");  
  
  FREE(data); FREE(data1); FREE(data2); 
  FREE(data_sphere);
  FREE(pack);
  
  MEM_STATS();
  
  return 0;
}

