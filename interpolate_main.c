#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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
#include "coarsen_surface.h"

#include "iso_modified.h"


#define SPHERE_DATA 1
#define DOUBLE_DATA 0
#define FLOAT_DATA  0

#define COARSEN_OUTPUT 0

#define MONITOR_MEMORY 0

#define VTK 0 /* vtk files output, in addition to Geomview (still working on it)*/



static void fill_from_double_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data);
      

int main (int argc, char * argv[])
{ 
   //if these options are changed recompile
  
  gboolean   color_plot = TRUE;
  gboolean   dbg_plot   = FALSE;
  gboolean   tetra      = FALSE; /* didn't prove to work well with
                                    interpolation */  				    
  gint       num_cells_to_peel = 3;
  
  GtsCartesianGrid g;
  gdouble isolevel, isolevel_fluid;
  
  GtsIsoCartesianFunc fill_function;
  
  gchar   base[256];
  gchar    fname1[256],  fnameg[256], fnamec[256];
  gdouble  *ddata_fluid, *ddata, *ddata2, *ddata_color;
  gint     n[3], nxyz;

  GtsSurface *surfg, *surf1, *surf1g, *surf12, *surf1_all, *surf_split1g, *surf_split12;
  GtsSurface *surfg2;
  GtsSurface  *surf2,*surf2g, *surf21, *surf2_all, *surf_split2g, *surf_split21;
  GNode      *tree_surfg;
  
  GSList     *triple_line;
  Pass_Data  *data;

  gint   i;
  //gint n1g, ng1;
  
  gchar  fname[256];
  FILE   *fp;
  guchar c;

  guint32 S0 = 0, S1 = 1, S2 = 2, SG = 4, S_DUMMY = 123456; //bitwise surface flags
  
  Data_Interpolate_Value  data_for_function, data_for_function_g;
 
 
   gdouble    d;
   gint   j, k, ind;
   gdouble  r_eff;
   gdouble  xc, yc, zc, h;
 
   gdouble  a, b, c1, x0, y0, z0, norm, dist, x, y, z;
   gdouble  center_color_value, color_value_min, color_value_max;
   
   gdouble  surf_area12, surf_area1, surf_area2, surf_areag, surf_area_core12;
   
   /* by default, coarsening will be done by number of edges */
   gboolean progressive = FALSE;
   gboolean log_cost = FALSE;
   StopOptions stop = COST; /* has worked better than NUMBER */
   CostOptions cost = COST_OPTIMIZED;
   MidvertexOptions mid = OPTIMIZED;
     
   
  /************************************************************************/
  
  /* if(argc < 4)
  {
      printf("\nUsage 1: interpolate fluid_data_file solid_data_file surf_output_fname ...");
      printf("\n      ...[color_data_file center_color_value]"); 
      printf("\nUsage 2: interpolate fluid_data_file solid_data_file surf_output_fname ...");
      printf("\n      ...[color_data_file min_color_value max_color_value]\n");
      
      return 1;
  } */
  
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
       
       if(argc == 4 )
           sprintf(base,"%s",argv[3]);    //basename for surface output (Geomview)
       else if (argc >= 5)
       {
           sprintf(fnamec,"%s",argv[3]);  //color info double array file
	   sprintf(base,"%s",argv[4]);    //basename for surface output (Geomview)
       }
  
       read_gdouble_array(&ddata_fluid,n,fname1,num_cells_to_peel);
       read_gdouble_array(&ddata,n,fnameg,num_cells_to_peel);
       
       if (argc >= 5) 
          read_gdouble_array(&ddata_color,n,fnamec,num_cells_to_peel);

       if (argc == 6)
       {
            center_color_value = atof(argv[5]);
	    printf("\ncenter_color_value %g", center_color_value); fflush(stdout);
       }
       else if(argc == 7)
       {
            color_value_min = atof(argv[5]);
	    color_value_max = atof(argv[6]);
       }
       
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
   else if( FLOAT_DATA )
  {
       /* level set method output assumed to be used:
        *  - fluid data is negative where the non-wetting fluid is,
	*  positive otherwise
	*  - mask(grain) data is negative in the pore space, positive in the
	*  grain space
       */
       sprintf(fname1,"%s",argv[1]);  //main fluid (fluid 1) float array file                         
       sprintf(fnameg,"%s",argv[2]);  //grain float array file
             
       if(argc == 4 )
           sprintf(base,"%s",argv[3]);    //basename for surface output (Geomview)
       else if (argc >= 5)
       {
           sprintf(fnamec,"%s",argv[3]);  //color info float array file
	   sprintf(base,"%s",argv[4]);    //basename for surface output (Geomview)
       }
 
         
       read_gfloat_array_into_gdouble(&ddata_fluid,n,fname1,num_cells_to_peel);
       read_gfloat_array_into_gdouble(&ddata,n,fnameg,num_cells_to_peel);
       
       if (argc >= 5) 
          read_gfloat_array_into_gdouble(&ddata_color,n,fnamec,num_cells_to_peel);
	  
       if (argc == 6)
       {
            center_color_value = atof(argv[5]);
	    printf("\ncenter_color_value %g", center_color_value); fflush(stdout);
       }
       else if (argc == 7)
       {
            color_value_min = atof(argv[5]);
	    color_value_max = atof(argv[6]);
       }
    
       
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
  
 
     // Define initial plane through (xc,yc,zc)
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
    
    sprintf(base,"interp_sphere");
  
 
  
     /**/ for(i = 0; i < nxyz; i++)
     {
       ddata[i] *= 0.02;
       ddata_fluid[i] *= 0.02;
       printf("\nhello");
       
       if( (ddata[i] > 0) && (ddata[i] > ddata_fluid[i]) ) ddata_fluid[i] = ddata[i];
     }  
     sprintf(base,"interp_hemisphere");   /**/     	
  }			    

      /* create opposite fluid phase */
   ddata2 = (gdouble *)MALLOC(nxyz * sizeof(gdouble));
  for(i = 0; i < nxyz; i++)
  {  
     gdouble tmp1 = -ddata_fluid[i];
     ddata2[i] = tmp1;
     //impose grain space masking
     if( (ddata[i] > ddata2[i]) ) ddata2[i] = ddata[i];
  }   

   /*******************************************************************************/
   
   /* Initialize surface structures */
   surf1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());

   surf2 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			    			       
   
   surfg = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			    			    			       			       
  
  /* fluid 1 surface */
  if( tetra )
    gts_isosurface_tetra(surf1, g, fill_function, ddata_fluid, isolevel_fluid);  
  else
    gts_isosurface_cartesian(surf1, g, fill_function, ddata_fluid, isolevel_fluid);  
    
  gts_surface_foreach_vertex(surf1,(GtsFunc)flag_initial,&S1);
 
  /* opposite fluid surface */
  if( tetra )
     gts_isosurface_tetra(surf2, g, fill_function, ddata2, isolevel_fluid);  
  else   
     gts_isosurface_cartesian(surf2, g, fill_function, ddata2, isolevel_fluid);
     
  gts_surface_foreach_vertex(surf2,(GtsFunc)flag_initial,&S2);  
    
  /* grain surface */
  if( tetra )
    gts_isosurface_tetra(surfg, g, fill_function, ddata, isolevel);
  else
    gts_isosurface_cartesian(surfg, g, fill_function, ddata, isolevel);
    
  gts_surface_foreach_vertex(surfg,(GtsFunc)flag_initial,&SG);
    
  
  /*interpolate ddata on surf1 */
  data_for_function.data = ddata_fluid;
  data_for_function.iso = isolevel_fluid;
  
  data_for_function.data_to_interpolate = ddata;
  data_for_function.iso_interpolate_data = isolevel;
  
  data_for_function.g = &g;
  data_for_function.flag_above = S1+SG;
  data_for_function.flag_below = S1+S2;
    
  gts_surface_foreach_vertex(surf1,(GtsFunc)gts_object_reset_reserved,NULL);
  
  if(argc >=5 )
  {
      data_for_function.data_to_interpolate_color = ddata_color;
      if (argc == 6)
      {
         data_for_function.color_min = 0.75*center_color_value;
	 data_for_function.color_max = 1.25*center_color_value;
	 data_for_function.use_color_minmax_input = 1;
	 printf("\nColor min and max set to 0.75 and 1.25 of the input value, %g, %g.",
	         data_for_function.color_min,data_for_function.color_max);
      }
      else if (argc == 7)
      {
         data_for_function.color_min = color_value_min;
	 data_for_function.color_max = color_value_max;
	 data_for_function.use_color_minmax_input = 1;
	 printf("\nColor min and max set input values, %g, %g.",
	         data_for_function.color_min,data_for_function.color_max);
      }
      else
      {
        data_for_function.color_min =  DBL_MAX;
        data_for_function.color_max = -DBL_MAX;
	data_for_function.use_color_minmax_input = 0;
	printf("\nColor min and max set from input data.");
      }
      gts_surface_foreach_vertex(surf1,(GtsFunc)vertex_interpolate_value_color,
                                                        &data_for_function);
							
      printf("\ncolor_min %g - GREEN, color_max %g - ORANGE",data_for_function.color_min,
                                                       data_for_function.color_max);							 
  }
  else
  {
      gts_surface_foreach_vertex(surf1,(GtsFunc)vertex_interpolate_value,
                                                        &data_for_function); 
  } 							
   /**interpolate ddata on surf2 **/
  data_for_function.data = ddata2;
  data_for_function.iso = isolevel_fluid;
  data_for_function.data_to_interpolate = ddata;
  data_for_function.iso_interpolate_data = isolevel;
  data_for_function.g = &g;
  data_for_function.flag_above = S2+SG;
  data_for_function.flag_below = S2+S1;
  
  
  gts_surface_foreach_vertex(surf2,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(surf2,(GtsFunc)vertex_interpolate_value,
                                                        &data_for_function); /**/
     					  
 
  /* Free input data */
  g_free(ddata_fluid);   g_free(ddata); g_free(ddata2); 
    
  gts_surface_foreach_face(surf1,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(surf2,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(surfg,(GtsFunc)flag_triangle_from_vert,NULL);
     
     surf_area1 = gts_surface_area(surf1);
     printf("\nsurf_area1  %g\n",surf_area1);
     surf_area2 = gts_surface_area(surf2);
     printf("\nsurf_area2  %g\n",surf_area2);
     surf_areag = gts_surface_area(surfg);
     printf("\nsurf_areag  %g\n",surf_areag);

/******/
  /* Get actual surfaces from all the flags*/
  surf1g =  get_surface_from_triangle_flags(surf1,S1+SG,S1);
  surf12 =  get_surface_from_triangle_flags(surf1,S1+S2,S1);
  surf1_all =  get_surface_from_triangle_flags(surf1,S1 + S2 + SG,S_DUMMY);
 
  surf_area_core12 = gts_surface_area(surf12);
  printf("\nsurf_area_core12  %g\n",surf_area_core12);
 
  //surfg2 = get_surface_from_triangle_flags(surfg,SG+S2,SG);
  
  split_surface_from_triangle_flags(surf1_all,S1+SG,S1+S2,isolevel_fluid,&surf_split1g,&surf_split12);
  
   /* Build bounding box trees for surfaces - used for efficient search*/
  tree_surfg = gts_bb_tree_surface (surfg);
  
   /* triple line edges */
  gts_surface_foreach_edge(surf_split12,(GtsFunc)flag_initial,&S0);
  gts_surface_foreach_edge(surf_split12,(GtsFunc)flag_edge_from_vertices,NULL);  

  triple_line = get_edgelist_from_edge_flags(surf_split12,S1 + S2 + SG,S_DUMMY);
  
  if (dbg_plot) gts_output_edgelist_oogl_spec(triple_line,&g,"interp_triple_line12.list");
  
  interpolate_measure_contact_angle_at_triple_line_edges(surf_split12,surfg,tree_surfg,triple_line,&g);
 
  
/******/
  //system("mv interp_testG1.list interp_testG1_1.list");
  //system("mv interp_test12.list interp_test12_1.list");
  surf2g =  get_surface_from_triangle_flags(surf2,S2+SG,S2);
  surf21 =  get_surface_from_triangle_flags(surf2,S2+S1,S2);
  surf2_all =  get_surface_from_triangle_flags(surf2,S1 + S2 + SG,S_DUMMY); 
  
  split_surface_from_triangle_flags(surf2_all,S2+SG,S2+S1,isolevel_fluid,&surf_split2g,&surf_split21);
    
  gts_surface_foreach_edge(surf_split21,(GtsFunc)flag_initial,&S0);
  gts_surface_foreach_edge(surf_split21,(GtsFunc)flag_edge_from_vertices,NULL);  

  triple_line = get_edgelist_from_edge_flags(surf_split21,S1 + S2 + SG,S_DUMMY);
  
  if (dbg_plot) gts_output_edgelist_oogl_spec(triple_line,&g,"triple_line21.list");
  
  interpolate_measure_contact_angle_at_triple_line_edges(surf_split21,surfg,tree_surfg,triple_line,&g);
								   
  if( color_plot )
  {
     /* Geomview files of (colored) surfaces are output by default.
      * Each file has local axes added so they can be superimposed at will.
      */ 

     if(argc >= 5)
     {
       /* output of interpolated values HAS to go first, note that 'gts_surface_write_oogl'
         resets the reserved field in vertices to NULL and thus voids interpolated values! */
       sprintf(fname,"%s_interp.list",base);
       fp = fopen(fname,"w");
       
       fprintf(fp,"LIST\n");
       print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);
       
       gts_surface_merge(surf12,surf_split12);       
       gts_surface_foreach_vertex(surf12,(GtsFunc)scale_interpolate_value_color,
                                                        &data_for_function);
       
       fprintf(fp,"\n{");
       gts_surface_write_oogl_color_vertices(surf12,fp);
       fprintf(fp,"\n}");	
	
       /* flush vertex color information */
       gts_surface_foreach_vertex(surf12,(GtsFunc)set_no_color,NULL);
       			
       gts_surface_merge(surf1g,surf_split1g);						
       Add_surf_only(surf1g,g,fp,"S1 & SG",set_no_color)
       fclose (fp);
     }

     gts_surface_foreach_face (surfg, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_grain.list",base);
     gts_output_surface_oogl_spec(surfg,&g,fname);

     if(VTK)
     {     
       sprintf(fname,"%s_grain.vtk",base);
       fp = fopen(fname,"w");
       gts_surface_write_vtk(surfg,fp);
       fclose(fp);
       
       sprintf(fname,"%s_fluid1.vtk",base);
       fp = fopen(fname,"w");
       gts_surface_write_vtk(surf1,fp);
       fclose(fp);
     }
     /* gts_surface_foreach_face (surfg2, (GtsFunc)set_green_color, NULL);
     sprintf(fname,"%s_g2.list",base);
     gts_output_surface_oogl_spec(surfg2,&g,fname); */
          
	  
     /* intersections file - for debugging
     sprintf(fname,"%s_intersect.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);

     //Add_surf(surf_split1g,g,fp,"S1 & SG",set_red_color,base,"_split_red.list")
     //Add_surf(surf_split12,g,fp,"SG & S2",set_no_color,base,"_split_grey.list")
     Add_surf(surf1g,g,fp,"S1 & SG",set_red_color,base,"_red.list")
     Add_surf(surf12,g,fp,"SG & S2",set_no_color,base,"_grey.list")
     Add_surf(surf1_all,g,fp,"SG & S1 & S2",set_yellow_color,base,"_yellow.list")
     */
     
     /* free surfaces not needed any more */
     gts_object_destroy (GTS_OBJECT (surf1));
     gts_object_destroy (GTS_OBJECT (surf2));
     gts_object_destroy (GTS_OBJECT (surfg));
     
     /* Fluid 1 surface, colored appropriately. */
     sprintf(fname,"%s_fluid1.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);     
     gts_surface_merge(surf1g,surf_split1g);
     gts_surface_merge(surf12,surf_split12);
     
     surf_area12 = gts_surface_area(surf12);
     printf("\nsurf_area12 %g\n",surf_area12);

     if (COARSEN_OUTPUT)
     
     {        
        gts_surface_coarsen_top(surf12,
	     progressive,cost,mid,log_cost,stop);
			 
	gts_surface_coarsen_top(surf1g, 
	     progressive,cost,mid,log_cost,stop);	 
     }
     
     Add_surf_only(surf12,g,fp,"S1 & S2",set_red_color)
     Add_surf_only(surf1g,g,fp,"S1 & SG",set_no_color)
     
     fclose (fp);
     
  
     if(VTK)
     {     
       sprintf(fname,"%s_mrg.vtk",base);
       fp = fopen(fname,"w");
       gts_surface_write_vtk(surf12,fp);
       fprintf(fp,"\n");
       gts_surface_write_vtk(surf1g,fp);
       fclose(fp);
     }
     
     /* free surfaces not needed any more */
     gts_object_destroy (GTS_OBJECT (surf12));
     gts_object_destroy (GTS_OBJECT (surf1g));
     gts_object_destroy (GTS_OBJECT (surf1_all));
     
     
     /* Fluid 2 surface, colored appropriately. */
     sprintf(fname,"%s_fluid2.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);     
     gts_surface_merge(surf2g,surf_split2g);
     gts_surface_merge(surf21,surf_split21);
     
     
    if (COARSEN_OUTPUT)
     {        
        gts_surface_coarsen_top(surf21,
	     progressive,cost,mid,log_cost,stop);
			 
	gts_surface_coarsen_top(surf2g, 
	     progressive,cost,mid,log_cost,stop);	 
     }
     
     Add_surf_only(surf21,g,fp,"S2 & S1",set_green_color)
     Add_surf_only(surf2g,g,fp,"S2 & SG",set_no_color)
     fclose (fp);     
  }     
  
  /* free surfaces not needed any more */
     gts_object_destroy (GTS_OBJECT (surf21));
     gts_object_destroy (GTS_OBJECT (surf2g));
     gts_object_destroy (GTS_OBJECT (surf2_all));
  
  if (MONITOR_MEMORY) MEM_STATS();
							 
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


