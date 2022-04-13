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
#include "segfl_util.h"
#include "viz_util.h"
#include "intersect_surf.h"
#include "cont_angle.h"
#include "mem.h"
#include "float.h"

#include "iso_modified.h"

#define SPHERE_DATA 0
#define DOUBLE_DATA 1
#define TOLERANCE 1e-10

#define INTERPOLATE 1 /* use basic interpolation along the edge where vertex is found */

void read_gdouble_array(gdouble **p_phi,gint *n,gchar *fname);
void write_gdouble_array(gdouble *phi,gint *n,gchar *fname);

static void fill_from_double_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data);
     

typedef struct _Data_Extrapolate_Fluid_Value{
    gdouble          *data_fluid;   //(packed) 3d array of fluid data
    gdouble          *data;         //(packed) 3d array of solid-grain data
    gdouble          iso;
    gdouble          iso_fluid;
    GtsCartesianGrid *g;
    guint32          flagG;           //flag for grain vertices
    guint32          flag1;          //flag for fluid 1 vertices
    guint32          flag2;          //flag for fluid 2 vertices
    GSList           *fluid_value;   //list of fluid values - not utilized for now
} Data_Extrapolate_Fluid_Value;


void vertex_extrapolate_fluid_value(GtsVertex *v, Data_Extrapolate_Fluid_Value *pdata)
{
    gint    i, j, k, i1, j1, k1, ind, test_x, test_y, test_z;
    gint    nx = (pdata->g)->nx, nxy = ((pdata->g)->nx) * ((pdata->g)->ny);    
    gint    shift_x  = 1, shift_y  =   nx, shift_z  =   nxy;
    gint    shift_2x = 2, shift_2y = 2*nx, shift_2z = 2*nxy;
    gint    shift_3x = 3, shift_3y = 3*nx, shift_3z = 3*nxy;
    
    gdouble   x,y,z;
    gdouble   test_sign_x, test_sign_y, test_sign_z;
    gdouble   *d, *f, d000, d100, d010, d001;
    gdouble   vplus1, vplus2, v1, v2;
    gdouble   fluid_value, *pfluid_value, test;
    gdouble   a, b, c, rho, one_minus_rho, f0, f1, f2;
    GtsCartesianGrid *g;
    guint32  flag;  
    
    g_return_if_fail (pdata->data_fluid != NULL);
    g_return_if_fail (pdata->data != NULL);
    
    /* shortcuts */
    d = pdata->data;
    f = pdata->data_fluid;
    g = pdata->g;
       
    x =  GTS_POINT(v)->x;
    y =  GTS_POINT(v)->y;
    z =  GTS_POINT(v)->z;
    //printf("\nHih (%g %g %g)",x,y,z);
    i =  floor(x - g->x);
    j =  floor(y - g->y);
    k = floor(z - g->z);
    
    i1 = ceil(x - g->x);
    j1 = ceil(y - g->y);
    k1 = ceil(z - g->z);
    
   
    if(i1 != i) test_x = 1; 
    else        test_x = 0;
    if(j1 != j) test_y = 1;
    else        test_y = 0;
    if(k1 != k) test_z = 1;
    else        test_z = 0;  

    ind = k*nxy + j*nx + i; //packed index in 3d array           
    d000 = d[ind] - pdata->iso;
    
    if ( i < (g->nx - 1) )  d100 = d[ind + shift_x] - pdata->iso;
    else                    d100 = d000;
    
    if ( j < (g->ny - 1) ) d010 = d[ind + shift_y] - pdata->iso;
    else                   d010 = d000;
     
    if ( k < (g->nz - 1) ) d001 = d[ind + shift_z] - pdata->iso;
    else                   d001 = d000;
    
    //printf("\ntests combined %d  d000 %g d100 %g d010 %g d001 %g",test_x+test_y+test_z,d000,d100,d010,d001);
    
    
    test_sign_x = d000*d100;
    test_sign_y = d000*d010;
    test_sign_z = d000*d001;
        
    //printf("\nHah!(%g %g %g) - (%d,%d,%d) - test_x %d %g test_y %d %g test_z %d %g",
    //                      x - g->x,y - g->y,z - g->z,i,j,k,
	//		  test_x,test_sign_x, test_y,test_sign_y, test_z,test_sign_z);
			    
    //if( test_x && (test_sign_x > 0)) printf(" sign_x");
    //if( test_y && (test_sign_y > 0)) printf(" sign_y");
    //if( test_z && (test_sign_z > 0)) printf(" sign_z");			    
    
    /* Determine on which cube edge is the vertex, and from which values to extrapolate */
    //if( ( d000 < 0. && d001 >= 0. ) || ( d000 >= 0. && d001 < 0. )  )
    if(test_z && (test_sign_z <= 0) )
    {
       v1 =  d000;    v2 = d001;
       if( v1 <= 0. ) //fluid vertex at position 000
       {
          f0 = f[ind];
	  
	  if (k >= 1) 
	  {
	     vplus1 = d[ind-shift_z] - pdata->iso;
	         f1 = f[ind-shift_z];
	  }  
	  else         
	     vplus1 = 0;
	     	
	  if (k >= 2)
	  {
	      vplus2 = d[ind-shift_2z] - pdata->iso;
	          f2 = f[ind-shift_2z]; 
	  }    
	  else 
	     vplus2 = 0;	  
       }
       else
       {
          f0 = f[ind+shift_z];
	  
	  if (k < g->nz-2)
	  {
	     vplus1 = d[ind+shift_2z] - pdata->iso;
	         f1 = f[ind+shift_2z];
	  }  
	  else         
	     vplus1 = 0;
	     
	  if (k < g->nz-3)
	  {
	      vplus2 = d[ind+shift_3z] - pdata->iso;
	          f2 = f[ind+shift_3z]; 
	  }    
	  else 
	     vplus2 = 0;		  
       }
       
       printf("\n%g %g %g -- k v1 %g v2 %g",x,y,z,v1,v2);         
    }
    else if (test_x && (test_sign_x <= 0) )
    //if( ( (d000 < 0.) && (d100 >= 0.) ) || ( (d000 >= 0.) && (d100 < 0.) ) )
    {    
       v1 =  d000;    v2 = d100;
       if( v1 <= 0. ) //fluid vertex at position 000
       {
          f0 = f[ind];  f_1  = f[ind+shift_x];
	  
	  if (i >= 1) 
	  {
	     vplus1 = d[ind-shift_x] - pdata->iso;
	         f1 = f[ind-shift_x];
	  }  
	  else         
	     vplus1 = 0;
	     	
	  if (i >= 2)
	  {
	      vplus2 = d[ind-shift_2x] - pdata->iso;
	          f2 = f[ind-shift_2x]; 
	  }    
	  else 
	     vplus2 = 0;	  
       }
       else
       {
          f0 = f[ind+shift_x]; f_1 = f[ind];
	  
	  if (i < g->nx-2)
	  {
	     vplus1 = d[ind+shift_2x] - pdata->iso;
	         f1 = f[ind+shift_2x];
	  }  
	  else         
	     vplus1 = 0;
	     
	  if (i < g->nx-3)
	  {
	      vplus2 = d[ind+shift_3x] - pdata->iso;
	          f2 = f[ind+shift_3x]; 
	  }    
	  else 
	     vplus2 = 0;		  
       } 
       printf("\n%g %g %g -- i=%d v1 %g v2 %g",x,y,z,i,v1,v2);     
    }
    else if (test_y&& (test_sign_y <= 0) ) 
    //if( ( d000 < 0. && d010>= 0. ) || ( d000 >= 0. && d010 < 0. ) )
    {
       v1 =  d000;    v2 = d010;
       if( v1 <= 0. ) //fluid vertex at position 000
       {
          f0 = f[ind]; f_1 =  f[ind+shift_y];
	  
	  if (j >= 1) 
	  {
	     vplus1 = d[ind-shift_y] - pdata->iso;
	         f1 = f[ind-shift_y];
	  }  
	  else         
	     vplus1 = 0;
	     	
	  if (j >= 2)
	  {
	      vplus2 = d[ind-shift_2y] - pdata->iso;
	          f2 = f[ind-shift_2y]; 
	  }    
	  else 
	     vplus2 = 0;	  
       }
       else
       {
          f0 = f[ind+shift_y]; f_1 = f[ind];
	  
	  if (j < g->ny-2)
	  {
	     vplus1 = d[ind+shift_2y] - pdata->iso;
	         f1 = f[ind+shift_2y];
	  }  
	  else         
	     vplus1 = 0;
	     
	  if (j < g->ny-3)
	  {
	      vplus2 = d[ind+shift_3y] - pdata->iso;
	          f2 = f[ind+shift_3y]; 
	  }    
	  else 
	     vplus2 = 0;		  
       }
       printf("\n%g %g %g -- j=%d v1 %g v2 %g",x,y,z,j,v1,v2);          
    }
    else
    {
      
      printf("\nBummer!(%g %g %g) - (%d,%d,%d) - test_x %d %g test_y %d %g test_z %d %g",
                          x - g->x,y - g->y,z - g->z,i,j,k,
			  test_x,test_sign_x, test_y,test_sign_y, test_z,test_sign_z);
      printf("\nx %g y %g z %g d000 %g d001 %g d010 %g d100 %g", x,y,z,d000,d001,d010,d100);			  
      
    }  
    
    /* Extrapolation */
    if( v1 <= 0. )  rho = v1/(v1 - v2);
    else            rho = v2/(v2 - v1);
    
    a = f0;
    
    pfluid_value = (gdouble *)malloc(sizeof(gdouble));
    
    
    //printf("\nv1 %g v2 %g vplus1 %g vplus2 %g",v1,v2,vplus1,vplus2);

/**/

    if ( INTERPOLATE ) fluid_value = f0 + rho*(f_1 - f0);
    else if ( vplus1 < 0. && vplus2 < 0. ) 
    {
	  b = -1.5*f0 + 2.0*f1 - 0.5*f2;
	  c =  0.5*f0 -     f1 + 0.5*f2; 
	  fluid_value = a - b * rho + c*rho*rho;
	  
	  if( fabs(fluid_value) < TOLERANCE) fluid_value = 0;
 	  
	  printf("\nf %g f0 %g f1 %g f2 %g", fluid_value,f0,f1,f2);
	  if(( ( fluid_value > f0) && (f0 < f1) ) || (( fluid_value > f0 ) && (f1 > f0))) printf(" test");
	  if( (( fluid_value > 0) && (f0 < 0)) || (( fluid_value < 0) && (f0 > 0)) )
	  {
	    printf(" test1");
	    fluid_value = 0;
	  }  
    }
    else 
    if ( vplus1 < 0.)
    {
	  b = f1 - a;
	  fluid_value = a - b * rho;
	  if( fabs(fluid_value) < TOLERANCE ) fluid_value = 0;
	  printf("\nf %g f0 %g f1 %g", fluid_value,f0,f1);
	  if(( ( fluid_value > f0) && (f0 < f1) ) || (( fluid_value > f0 ) && (f1 > f0) ) ) printf(" test");
          if( (( fluid_value > 0) && (f0 < 0)) || (( fluid_value < 0) && (f0 > 0)) ) 
	  {
	    printf(" test1");
	    fluid_value = 0;
	  }  	 
    } 
    else   
    {
          fluid_value = a;
	  printf("\nf %g f0 %g", fluid_value,f0);
    }
            
    test = fluid_value - pdata->iso_fluid;
    if ( test > 0)
      flag = pdata->flag2 + pdata->flagG;
    else //test is effectively 0 or it's negative which points to fluid 1
      flag = pdata->flag1 + pdata->flagG;
      
    printf("  flag %d",flag);  
    GTS_OBJECT(v)->flags = flag;
    *pfluid_value = fluid_value;
    GTS_OBJECT(v)->reserved = pfluid_value;        
}


int main (int argc, char * argv[])
{ 
   //if these options are changed recompile
  gboolean   verbose =    FALSE;
  
  gboolean   color_plot = TRUE;
  gboolean   dbg_plot =   FALSE;
  gboolean   web_plot =   FALSE;
  gboolean   print_pointwise = FALSE;  
  
  GtsCartesianGrid g;
  gdouble isolevel, isolevel_fluid;
  
  GtsIsoCartesianFunc fill_function;
  
  gchar   base[256];
  gchar    fname1[256],  fnameg[256];
  gdouble  *ddata_fluid, *ddata;
  gint     n[3], nxyz;
  

  GtsSurface *surfg, *surfg1, *surfg2, *surfg_all, *surf_split1, *surf_split2, *surf1;
  //GNode      *tree1, *tree2, *treeg;
  //GtsSurface *surf12,  *surfg1, *surfg2, *surf1g, *surf21, *surf2g;
  //GtsSurface *surf_1only, *surf_2only, *surf_gonly;

  gint   i, blob_id = 0;
  //gint n1g, ng1;
  
  gchar  fname[256];
  FILE   *fp;
  guchar c;

  guint32 S0 = 0, S1 = 1, S2 = 2, SG = 4, S_DUMMY = 123456; //bitwise surface flags
  
  Data_Extrapolate_Fluid_Value  data_for_function;
 
 
   gdouble    d;
   gint   j, k, ind;
   gdouble  r_eff;
   gdouble  xc, yc, zc, h;
 
   gdouble  a, b, c1, x0, y0, z0, norm, dist, x, y, z;
 
  /************************************************************************/
  
 MEM_START();
 
  /* Read in double data*/
  if( DOUBLE_DATA )
  {
       /* level set method output assumed to be used:
        *  - fluid data is negative where the non-wetting fluid is,
	*  positive otherwise
	*  - mask(grain) data is negative in the pore space, positive in the
	*  grain space
       */
       sprintf(fname1,"%s",argv[1]);  //main fluid (fluid 1) double array file                         
       sprintf(fnameg,"%s",argv[2]);  //grain double array file
       sprintf(base,"%s",argv[3]);    //basename for surface output (Geomview)
       
       if (web_plot || print_pointwise ) blob_id     = atoi(argv[4]);
  
       read_gdouble_array(&ddata_fluid,n,fname1);
       read_gdouble_array(&ddata,n,fnameg); 
       
       nxyz = n[0]*n[1]*n[2];              
	
       /* this creates segmented type of data, if desired *
       for(i = 0; i < nxyz; i++)
       {          
	   if(ddata[i]  > 0)  
	        ddata[i] =  0.5;
	   else ddata[i] = -0.5;
	   
	   if(ddata_fluid[i] > 0) 
	        ddata_fluid[i] =  0.5; 
	   else ddata_fluid[i] = -0.5;	       
       }	
       ***/
	
       isolevel = 0.0;
       isolevel_fluid = 0.0;
       fill_function = fill_from_double_data;
       
          /* Initialize grid */
       g.dx = g.dy = g.dz = 1.0; //voxel length assumed 1.0 for now
       g.x = g.y = g.z = -0.5;   //center of the first voxel assumed at (0,0,0)
                             //hence volume boundary starts at (-.5,-.5,-.5)
       g.nx = n[0]; g.ny = n[1]; g.nz = n[2];

   }    
   else /* SPHERE */
   {   
       isolevel = 0.0;
       isolevel_fluid = 0.0;
       fill_function = fill_from_double_data;
       
       sprintf(base,"sphere");
       
      /*
      * create double data_sphere with diameter 'd'
      * fluid data is negative at voxels inside the sphere, otherwise positive
     */
      d = 7.0;
      
       /* Initialize grid */
      g.dx = g.dy = g.dz = 1.0; //voxel length assumed 1.0 for now
      g.x = g.y = g.z = -0.5;   //center of the first voxel assumed at (0,0,0)
                            //hence volume boundary starts at (-.5,-.5,-.5)
      g.nx = g.ny = g.nz = d+6;
  
      nxyz = g.nx * g.ny * g.nz;
      n[0] = g.nx;
      n[1] = g.ny;
      n[2] = g.nz;
      ddata =       (gdouble *)MALLOC(nxyz * sizeof(gdouble));
      ddata_fluid = (gdouble *)MALLOC(nxyz * sizeof(gdouble));
  
     
      xc = yc = zc = 0.5*d + 2.2; //sphere center
  
      r_eff = 0.5*d; //sphere radius
      ind = 0;
      for(k = 0; k < g.nz; k++) 
       for(j = 0; j < g.ny; j++)
        for(i = 0; i < g.nx; i++, ind++)
        {	
	   x = g.x + i*g.dx; 
	   y = g.y + j*g.dy; 
	   z = g.z + k*g.dz;  
	  dist = (x - xc)*(x - xc) + (y - yc)*(y - yc) + (z - zc)*(z - zc);
	  dist = sqrt(dist);
	  
          ddata_fluid[ind] = dist - r_eff;
        }
  
 
     // Define intial plane through (xc,yc,zc)
     //default value
     a = 0; b = 1; c1 = 1;  //(a,b,c) is plane normal vector
  
  
     norm = sqrt(a*a + b*b + c1*c1);
     //normalize
     a = a/norm; b = b/norm; c1 = c1/norm;
  
//plane moves if normal direction away from sphere center (xc,yc,zc)
//    for(h = 0; h < r; h++)
//    {
//       printf("\nh %g",h); fflush(stdout);
     
       //new plane position
       h = 0;
       
       x0 = xc + h*a;
       y0 = yc + h*b;
       z0 = zc + h*c1;
      
      ind = 0;
      for(k = 0; k < g.nz; k++) 
       for(j = 0; j < g.ny; j++)
        for(i = 0; i < g.nx; i++, ind++)
        {
           x = g.x + i*g.dx; 
	   y = g.y + j*g.dy; 
	   z = g.z + k*g.dz;  
	  ddata[ind] = a*(x - x0) + b*(y - y0) + c1*(z - z0);
        }
     
//   }
    }
       

   /* Initialize (grain) surface structure */
   surfg = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
			    

   surf1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			    
			    
   			    

  /* 
  * Produce triangulated surfaces by marching cubes algorithm.
  * Note that each triangle is oriented so that its normal by default points
  * to larger data values (i.e. 1 in our case).
  * Therefore surf1 normals point outwards from fluid 1 (surf2, surfg similarly)
  * Modification extrapolates fluid data onto the grain surface vertices:
  * each vertex should be flagged to belong to fluid 1 (S1+SG) or fluid2 (S2+SG)
  * 
   gts_isosurface_cartesian_modified(surfg, g, fill_function, ddata,  
                                     ddata_fluid, isolevel, isolevel_fluid); */
  
  
  printf("\nIsosurfacing"); fflush(stdout);
  gts_isosurface_cartesian_copy(surfg, g, fill_function, ddata, isolevel);
  
  gts_surface_foreach_vertex(surfg,(GtsFunc)flag_initial,&SG);
  
  /* this is to get a half-sphere as a fluid blob */
  for(i=0; i < nxyz; i++)
  {   
     if( (ddata[i] > 2.0) && (ddata[i] > ddata_fluid[i]) ) 
        ddata_fluid[i] = ddata[i];
  }
  
  gts_isosurface_cartesian(surf1, g, fill_function, ddata_fluid, isolevel_fluid);
  
  
  data_for_function.data = ddata;
  data_for_function.iso = isolevel;
  data_for_function.data_fluid = ddata_fluid;
  data_for_function.iso_fluid = isolevel_fluid;
  data_for_function.g = &g;
  data_for_function.flag1 = S1;
  data_for_function.flag2 = S2;
  data_for_function.flagG = SG;
  
  printf("\nExtrapolation"); fflush(stdout);
  gts_surface_foreach_vertex(surfg,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(surfg,(GtsFunc)vertex_extrapolate_fluid_value,
                                   &data_for_function);	
  
  
  /* temporary output for contact_angle 
  for(i=0; i < nxyz; i++)
  {   if( ddata[i] > ddata_fluid[i]) ddata_fluid[i] = ddata[i]; } 
     
  write_gdouble_array(ddata_fluid,n,"sphere_fluid_data");
  write_gdouble_array(ddata,n,"sphere_mask_data");     */

  /* Free input data */
  g_free(ddata_fluid);   g_free(ddata);  
  
   
  /* Display summary information about surfaces */
//  if (verbose)
  {
    fprintf(stderr,"\n#GRAIN surface SG");
    gts_surface_print_stats (surfg, stderr);
  }
  fflush(stderr);

  gts_surface_foreach_face(surfg,(GtsFunc)flag_triangle_from_vert,NULL);

  /* Get actual surfaces from all the flags*/
  surfg1 =  get_surface_from_triangle_flags(surfg,SG+S1,SG);
  surfg2 =  get_surface_from_triangle_flags(surfg,SG+S2,SG);
  surfg_all =  get_surface_from_triangle_flags(surfg,S1 + S2 + SG,S_DUMMY);
 
  split_surface_from_triangle_flags(surfg_all,SG+S1,SG+S2,isolevel_fluid,&surf_split1,&surf_split2);	
								   
  if( color_plot )
  {
     /* Geomview files of (colored) surfaces are output by default.
      * Each file has local axes added so they can be superimposed at will.
      */ 

     gts_surface_foreach_face (surfg, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_grain.list",base);
     gts_output_surface_oogl_spec(surfg,&g,fname);
     
     gts_surface_foreach_face (surf1, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_fluid.list",base);
     gts_output_surface_oogl_spec(surf1,&g,fname);

     //intersections file
     sprintf(fname,"%s_intersect.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);

     Add_surf(surf_split1,g,fp,"S1 & SG",set_magenta_color,base,"_split_magenta.list")
     Add_surf(surf_split2,g,fp,"SG & S2",set_cyan_color,base,"_split_cyan.list")
     Add_surf(surfg1,g,fp,"S1 & SG",set_magenta_color,base,"_magenta.list")
     Add_surf(surfg2,g,fp,"SG & S2",set_cyan_color,base,"_cyan.list")
     Add_surf(surfg_all,g,fp,"SG & S1 & S2",set_yellow_color,base,"_yellow.list")
     
     fclose (fp);
  }
  
  
  
  
  MEM_STATS();
							 
  return 0;
}



void read_gdouble_array(gdouble **p_phi,gint *n,gchar *fname)
{
   FILE *fp;
   gint  nxyz;
   gdouble *phi;
   
   fp = fopen(fname,"r");
   fread(n,sizeof(gint),3,fp); //read dimensions n[0],n[1],n[2]
   
   nxyz = n[0]*n[1]*n[2];
   phi = (gdouble *)g_malloc(nxyz*sizeof(gdouble));
   fread(phi,sizeof(gdouble),nxyz,fp);
   fclose(fp);
   
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


static void fill_from_double_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data)
{
  guint	   i, j;
  gdouble  *pdata;
  gint     nxy = (g.nx)*(g.ny);
  
  //printf("\nk %d nxy %d",k, nxy);
  pdata = (gdouble *)data + k*nxy;
  for( j = 0; j < g.ny; j++ )
  {
	for( i = 0; i < g.nx; i++, pdata++ )
	{
		f[i][j] = *pdata;
        }
  }
}			
