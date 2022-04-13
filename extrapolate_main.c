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
#include "double_data_util.h"
#include "iso_modified.h"

#define SPHERE_DATA 0
#define DOUBLE_DATA 0
#define FLOAT_DATA 1


static void fill_from_double_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data);
     


int main (int argc, char * argv[])
{ 
   //if these options are changed recompile  
  gboolean   color_plot = TRUE;
  gint       num_cells_to_peel = 0;
  
  GtsCartesianGrid g;
  gdouble isolevel, isolevel_fluid;
  
  GtsIsoCartesianFunc fill_function;
  
  gchar   base[256];
  gchar    fname1[256],  fnameg[256];
  gdouble  *ddata_fluid, *ddata;
  gint     n[3], nxyz;
  

  GtsSurface *surfg, *surfg1, *surfg2, *surfg_all, *surf_splitG1, *surf_splitG2, *surf1;
  
  GtsSurface *surf12;

  gint   i;
  //gint n1g, ng1;
  
  gchar  fname[256];
  FILE   *fp;
  guchar c;

  guint32 S0 = 0, S1 = 1, S2 = 2, SG = 4, S_DUMMY = 123456; //bitwise surface flags
  
  Data_Extrapolate_Fluid_Value  data_for_function;
  Data_Interpolate_Value        data_for_function_interp;
 
   gdouble    d;
   gint   j, k, ind;
   gdouble  r_eff;
   gdouble  xc, yc, zc, h;
 
   gdouble  a, b, c1, x0, y0, z0, norm, dist, x, y, z;
 
 
  GSList     *triple_line;
  /************************************************************************/
  
 MEM_START();
 
  /* Read in double data*/
  if( FLOAT_DATA )
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
         
       read_gfloat_array_into_gdouble(&ddata_fluid,n,fname1,num_cells_to_peel);
       read_gfloat_array_into_gdouble(&ddata,n,fnameg,num_cells_to_peel);

       
       nxyz = n[0]*n[1]*n[2];              
	for(i = 0; i < nxyz; i++) ddata[i] = -ddata[i];
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
  else if( DOUBLE_DATA )
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
         
       read_gdouble_array(&ddata_fluid,n,fname1,num_cells_to_peel);
       read_gdouble_array(&ddata,n,fnameg,num_cells_to_peel);      
       
       
       nxyz = n[0]*n[1]*n[2];              
	for(i = 0; i < nxyz; i++) ddata[i] = -ddata[i];
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

     sprintf(base,"extrapol_sphere");
  
     
    /***/ for(i = 0; i < nxyz; i++) 
       if( (ddata[i] > ddata_fluid[i]) ) ddata_fluid[i] = ddata[i];
     sprintf(base,"extrapol_hemisphere");
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
  gts_isosurface_cartesian(surfg, g, fill_function, ddata, isolevel);  
  gts_surface_foreach_vertex(surfg,(GtsFunc)flag_initial,&SG);
   
  gts_isosurface_cartesian(surf1, g, fill_function, ddata_fluid, isolevel_fluid);
  gts_surface_foreach_vertex(surf1,(GtsFunc)flag_initial,&S1);
  
  
  /* Extrapolate fluid data onto grain surface */
  data_for_function.data = ddata;
  data_for_function.iso = isolevel;
  data_for_function.data_fluid = ddata_fluid;
  data_for_function.iso_fluid = isolevel_fluid;
  data_for_function.g = &g;
  data_for_function.flag1 = S1;
  data_for_function.flag2 = S2;
  data_for_function.flagG = SG;
    
  gts_surface_foreach_vertex(surfg,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(surfg,(GtsFunc)vertex_extrapolate_fluid_value,
                                   &data_for_function);	
  
  
  /* Interpolate ddata on surf1 */
  data_for_function_interp.data = ddata_fluid;
  data_for_function_interp.iso = isolevel_fluid;
  
  data_for_function_interp.data_to_interpolate = ddata;
  data_for_function_interp.iso_interpolate_data = isolevel;
  
  data_for_function_interp.g = &g;
  data_for_function_interp.flag_above = S1+SG;
  data_for_function_interp.flag_below = S1+S2;
    
  gts_surface_foreach_vertex(surf1,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(surf1,(GtsFunc)vertex_interpolate_value,
                                                        &data_for_function_interp); 
  
  /* temporary output for contact_angle 
  for(i=0; i < nxyz; i++)
  {   if( ddata[i] > ddata_fluid[i]) ddata_fluid[i] = ddata[i]; } 
     
  write_gdouble_array(ddata_fluid,n,"sphere_fluid_data");
  write_gdouble_array(ddata,n,"sphere_mask_data");     */

  /* Free input data */
  g_free(ddata_fluid);   g_free(ddata);  
  
  gts_surface_foreach_face(surfg,(GtsFunc)flag_triangle_from_vert,NULL);

  /* Get actual surfaces from all the flags*/
  surfg1 =  get_surface_from_triangle_flags(surfg,SG+S1,SG);
  surfg2 =  get_surface_from_triangle_flags(surfg,SG+S2,SG);
  surfg_all =  get_surface_from_triangle_flags(surfg,S1 + S2 + SG,S_DUMMY);
 
  split_surface_from_triangle_flags(surfg_all,SG+S1,SG+S2,isolevel_fluid,&surf_splitG1,&surf_splitG2);
  
  gts_surface_foreach_edge(surf_splitG1,(GtsFunc)flag_initial,&S0);
  gts_surface_foreach_edge(surf_splitG1,(GtsFunc)flag_edge_from_vertices,NULL);  

  triple_line = get_edgelist_from_edge_flags(surf_splitG1,S1 + S2 + SG,S_DUMMY);
  
  gts_output_edgelist_oogl_spec(triple_line,&g,"triple_line.list");
  
   
  gts_surface_foreach_face(surf1,(GtsFunc)flag_triangle_from_vert,NULL);
  surf12 =  get_surface_from_triangle_flags(surf1,S1+S2,S1);
  
  extrapolate_measure_contact_angle_at_triple_line_edges(surf_splitG1,surf12,triple_line,&g);
 
    								   
  if( color_plot )
  {
     /* Geomview files of (colored) surfaces are output by default.
      * Each file has local axes added so they can be superimposed at will.
      */ 

     gts_surface_foreach_face (surf1, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_fluid1.list",base);
     gts_output_surface_oogl_spec(surf1,&g,fname);
     
     gts_surface_foreach_face (surf12, (GtsFunc)set_red_color, NULL);
     sprintf(fname,"%s_f12.list",base);
     gts_output_surface_oogl_spec(surf12,&g,fname);

     //intersections file
     sprintf(fname,"%s_intersect.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);

     //Add_surf(surf_splitG1,g,fp,"S1 & SG",set_magenta_color,base,"_split_magenta.list")
     //Add_surf(surf_split2,g,fp,"SG & S2",set_cyan_color,base,"_split_cyan.list")
     Add_surf(surfg1,g,fp,"S1 & SG",set_magenta_color,base,"_magenta.list")
     Add_surf(surfg2,g,fp,"SG & S2",set_cyan_color,base,"_cyan.list")
     Add_surf(surfg_all,g,fp,"SG & S1 & S2",set_yellow_color,base,"_yellow.list")
     
     sprintf(fname,"%s_grain.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);     
     gts_surface_merge(surfg1,surf_splitG1);
     gts_surface_merge(surfg2,surf_splitG2);
     Add_surf(surfg1,g,fp,"S1 & SG",set_magenta_color,base,"_magenta.list")
     Add_surf(surfg2,g,fp,"SG & S2",set_cyan_color,base,"_cyan.list")
     
     fclose (fp);
  }
  
  
  
  
  MEM_STATS();
							 
  return 0;
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
