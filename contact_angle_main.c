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

#include "double_data_util.h"

static void fill_from_data(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data);

static void fill_from_double_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data);
      
void  create_test_example(gdouble **pdata1, gdouble **pdata2,gdouble **pdata,Cube_Info *);
      

/* IF CODE IS COMPILED WITH seg_data == TRUE 
*  Run as './contact_angle segfl_name1 segfl_name2 base'
*  or    './contact_angle segfl_name1 segfl_name2 base blob_id' if
*         the case that web_plot = TRUE

*  If test is set to TRUE and recompile, a simple test example executable
*  './contact_angle' is created. 
*
*
*  ELSE IF INPUT DATA ARE DOUBLE ARRAYS (COMPILE CODE WITH seg_data == 0)
*   Run as './contact_angle fname1 fname_mask base'
*   or    './contact_angle  fname1 fname_mask base blob_id' if
*                      (fluid1 fluid2 grain)
*   Note that fluid2 data will be created from fluid1 and mask data.		      
*		      
*  web_plot = TRUE

*  
*/


/* 2022/04 change
   Usage options have a number of modes.
   
   contact_angle 
   contact_angle 1   
        run a test example (simple geometry with three "block" shaped regions for fluid1, 2 and grain
   
   contact_angle 2 fname1 fname_mask base
        run a code with float data arrays with level set functions describing data
	fname1     = level set function where data is < 0 where fluid 1 is located
	fname_mask = level set function where data is < 0 where pore space is located
	base       = base file name for Geomview visualization results
	blob_id    = optional input
	Note that fluid 2 is where we have pore space and no fluid 1.
   
   contact_angle 3 fname1 fname_mask base
        run a code with double data arrays with level set functions describing data
	fname1     = level set function where data is < 0 where fluid 1 is located
	fname_mask = level set function where data is < 0 where pore space is located
	base       = base file name for Geomview visualization results
	blob_id    = optional input
	Note that fluid 2 is where we have pore space and no fluid 1.
  
   	
   contact_angle 4 segfl_name1 segfl_name2 base
   contact_angle 4 segfl_name1 segfl_name2 blob_id
   
        segfl_name1, segfl_name2 are segmented files from 3DMA-Rock	
	base = base file name for Geomview visualization results
	blob_id  = if provided, it assumes that
	
   TO BE IMPLEMENTED PROPERLY	
   contact_angle 5 NX NY NZ segfl_uchar base
       NX, NY, NZ positive integers denoting the 3D array dimensions
       segfl_uchar is assumed a binary array with unsigned character data that is set to 
                   0,1,2 where grain, fluid1, fluid2 phases are respectively
*/


int main (int argc, char * argv[])
{ 
  // the following options are set based on the input arguments 
  gboolean   test =             FALSE;
  gboolean   seg_data_3dma =    FALSE;
  gboolean   seg_data_uchar=    FALSE; 
  gboolean   double_data =      FALSE;
  gboolean   float_data =       FALSE;
  gboolean   web_plot =         FALSE;
  
  
  // To change the following options, you need to recompile the code
  gboolean   color_plot = TRUE;
  gboolean   dbg_plot =   FALSE;
  gboolean   print_pointwise = FALSE;
  gboolean   verbose =    FALSE;
  gboolean   tetra =      FALSE;
   
  guint      num_cells_to_peel = 0; 
  
  GtsCartesianGrid g;
  gdouble isolevel;
  
  GtsIsoCartesianFunc fill_function;
  
  gchar   segfl_name1[256], segfl_name2[256], base[256], segfl_name[256];
  guchar  *data1, *data2, *data, c;
  Cube_Info   cube;
  
  gchar    fname1[256],  fnameg[256];
  gdouble  *ddata1, *ddata2, *ddata, tmp1, tmp2;
  gint     n[3], nxyz;
  

  GtsSurface *surf1, *surf2, *surfg;
  GNode      *tree1, *tree2, *treeg;
  GtsSurface *surf12,  *surfg1, *surfg2, *surf1g, *surf21, *surf2g;
  GtsSurface *surf_1only, *surf_2only, *surf_gonly;

  gint   i, mode, blob_id = 0;
  //gint n1g, ng1;
  
  gchar  fname[256];
  FILE   *fp;

  guint32 all, S1 = 1, S2 = 2, SG = 4, S_DUMMY = 0; //bitwise surface flags
  
  gdouble  area_12, area_g1, area_g2, area_1only, area_2only, area_gonly;
  gdouble  area_1, vol_1;

  Pass_Data_vca  vca1, vca2;
  GSList   *master_vert_list;
  
  /************************************************************************/
  
  //MEM_START();
 
  //printf("\n %d %d %d %d %d",test, float_data, double_data, seg_data,web_plot);
  //fflush(stdout);
  
  if (argc <= 1)  test = TRUE;
  else
  {
      mode = atoi(argv[1]);
      printf("\nMode %d",mode);
      fflush(stdout);
      
      switch (mode)
      {
         case 1:
      	          test = TRUE;
		  break;
	 case 2:    
                  float_data = TRUE;
		  break;
         case 3:
	 	  double_data = TRUE;
		  break;
	 case 4:
	 	  seg_data_3dma = TRUE;
		  if (argc == 6) web_plot = TRUE; //if blob_id provided, set plotting accordingly
		  break;
	
	 case 5:  seg_data_uchar = TRUE;
	          break; 	  
	 default:
	 	  test = TRUE; 
		  break;
		  		  			    
      }   
  }
  
  // printf("\n %d %d %d %d %d",test, float_data, double_data, seg_data,web_plot);
  
  /* Read in data*/
  if( test )
  {
      printf("\ntest mode, creating data.");
      create_test_example(&ddata1,&ddata2,&ddata,&cube);
      
      sprintf(base,"test");
      n[0] = cube.nx;
      n[1] = cube.ny;
      n[2] = cube.nz;
      isolevel = 0.5;
      fill_function = fill_from_double_data;
  }
  else if( seg_data_3dma )
  {
      printf("\nseg_data mode, 3DMA-Rock segmented data assumed as input.");
      sprintf(segfl_name1,"%s",argv[2]); //main fluid (fluid 1) segmented file
      sprintf(segfl_name2,"%s",argv[3]); //second fluid (fluid 2) segmented file
      sprintf(base,"%s",argv[4]);        //basename for surface output (Geomview)
      if (web_plot || print_pointwise )  blob_id     = atoi(argv[5]);

      /*
      * Read segmented file information.
      * Fluid segmented files assumed to have 0 where fluid is, 1 otherwise.
      * Structure 'cube' is used for historical reasons (3DMA)
      * TO DO: initialize grid directly from reading seg files.
      */
      
      printf("\n%s %s %s %d",segfl_name1, segfl_name2,base,blob_id);
      fflush(stdout);
      read_segfl_top(segfl_name1,&data1,&cube);
      read_segfl_top(segfl_name2,&data2,&cube);
  

     /*
     * Create 'data' - 0 where grain is, 1 otherwise.
     */
     data = (guchar *)g_malloc(cube.nxyz*UCSZ );
     for(i = 0; i < cube.nxyz; i++)
     {
	 if( (data1[i] == 1) && (data2[i] == 1) ) data[i] = 0;
	 else                                     data[i] = 1;
     }

     isolevel = 0.5;
     fill_function = fill_from_data;
     n[0] = cube.nx;
     n[1] = cube.ny;
     n[2] = cube.nz;
  }
  else if( seg_data_uchar )
  {
      printf("\nseg_data_uchar mode");
      n[0] = atoi(argv[2]);
      n[1] = atoi(argv[3]);
      n[2] = atoi(argv[4]);
      
      sprintf(segfl_name,"%s",argv[5]); 
      sprintf(base,"%s",argv[6]);        //basename for surface output (Geomview)
      if (web_plot || print_pointwise )  blob_id     = atoi(argv[7]);

      /*
      * Read segmented file information.
      * Fluid segmented files assumed to have 0 where fluid is, 1 otherwise.
      * Structure 'cube' is used for historical reasons (3DMA)
      * TO DO: initialize grid directly from reading seg files.
      */
      
      printf("\n%d %d %d %s %s %s %d",n[0],n[1],n[2],segfl_name,base,blob_id);
      fflush(stdout);
 
     /*
     * Create 'data' - 0 where grain is, 1 otherwise.
     */
     //data = (guchar *)g_malloc(cube.nxyz*UCSZ );
     //for(i = 0; i < cube.nxyz; i++)
     //{
     //	 if( (data1[i] == 1) && (data2[i] == 1) ) data[i] = 0;
     //	 else                                     data[i] = 1;
     //}

     isolevel = 0.5;
     fill_function = fill_from_data;
  }
  
  else if(double_data)
  {
       printf("\ndouble data mode, LSMPQS 'double' level set function input data assumed");
       sprintf(fname1,"%s",argv[2]);  //main fluid (fluid 1) double array file
       sprintf(fnameg,"%s",argv[3]);  //grain double array file
       sprintf(base,"%s",argv[4]);  //basename for surface output (Geomview)
       if (web_plot || print_pointwise ) blob_id     = atoi(argv[5]);
  
       read_gdouble_array(&ddata1,n,fname1,num_cells_to_peel);
       read_gdouble_array(&ddata,n,fnameg,num_cells_to_peel);
       
       /* create fluid2 data from fluid1 and grain */
       nxyz = n[0]*n[1]*n[2];
       ddata2 = (gdouble *)g_malloc(nxyz*sizeof(gdouble));
       
       for(i = 0; i < nxyz; i++)
       {
           tmp1 = -ddata1[i];
	   tmp2 = ddata[i];
	    
	   ddata2[i] = ( tmp1 > tmp2 ) ? tmp1 : tmp2;
	   
	   //if(ddata[i]  > 0)  ddata[i] = 0.5; else ddata[i]  = -0.5;
	   //if(ddata1[i] > 0) ddata1[i] = 0.5; else ddata1[i] = -0.5;
	   //if(ddata2[i] > 0) ddata2[i] = 0.5; else ddata2[i] = -0.5;	    
       }
       
       //lsm double data is positive where grain space is, need to negate
       for(i = 0; i < nxyz; i++) ddata[i] = -ddata[i];
       
       isolevel = 0.0;
       fill_function = fill_from_double_data;
   }
   else
   {
       printf("\nfloat data mode, LSMPQS 'float' level set function input data assumed");
       sprintf(fname1,"%s",argv[2]);  //main fluid (fluid 1) float array file
       sprintf(fnameg,"%s",argv[3]);  //grain float array file
       sprintf(base,"%s",argv[4]);  //basename for surface output (Geomview)
       if (web_plot || print_pointwise ) blob_id     = atoi(argv[5]);
       
       read_gfloat_array_into_gdouble(&ddata1,n,fname1,num_cells_to_peel);
       read_gfloat_array_into_gdouble(&ddata,n,fnameg,num_cells_to_peel);
       
              /* create fluid2 data from fluid1 and grain */
       nxyz = n[0]*n[1]*n[2];
       ddata2 = (gdouble *)g_malloc(nxyz*sizeof(gdouble));
       
       for(i = 0; i < nxyz; i++)
       {
           tmp1 = -ddata1[i];
	   tmp2 = ddata[i];
	    
	   ddata2[i] = ( tmp1 > tmp2 ) ? tmp1 : tmp2;
	   
	   //if(ddata[i]  > 0)  ddata[i] = 0.5; else ddata[i]  = -0.5;
	   //if(ddata1[i] > 0) ddata1[i] = 0.5; else ddata1[i] = -0.5;
	   //if(ddata2[i] > 0) ddata2[i] = 0.5; else ddata2[i] = -0.5;	    
       }
       
       //LSMPQS data is positive where grain space is, need to negate
       for(i = 0; i < nxyz; i++) ddata[i] = -ddata[i];
       
       isolevel = 0.0;
       fill_function = fill_from_double_data;
 
   }
      
       
   /* Initialize grid */
   g.dx = g.dy = g.dz = 1.0; //voxel length assumed 1.0 for now
   g.x = g.y = g.z = -0.5;   //center of the first voxel assumed at (0,0,0)
                             //hence volume boundary starts at (-.5,-.5,-.5)
   g.nx = n[0]; g.ny = n[1]; g.nz = n[2];

   /* Initialize surface structures */
   surf1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
   surf2 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
   surfg = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());

  /* 
  * Produce triangulated surfaces by marching cubes algorithm.
  * Note that each triangle is oriented so that its normal by default points
  * to larger data values (i.e. 1 in our case).
  * Therefore surf1 normals point outwards from fluid 1 (surf2, surfg similarly). 
  */

  if(seg_data_3dma)
  {
     if (tetra)
     {
      gts_isosurface_tetra (surf1, g, fill_function, data1, isolevel);
      gts_isosurface_tetra (surf2, g, fill_function, data2, isolevel);
      gts_isosurface_tetra (surfg, g, fill_function, data,  isolevel);
     }
     else
     {
      gts_isosurface_cartesian (surf1, g, fill_function, data1, isolevel);
      gts_isosurface_cartesian (surf2, g, fill_function, data2, isolevel);
      gts_isosurface_cartesian (surfg, g, fill_function, data,  isolevel);
     }

     /* Free data */
     g_free(data1);    g_free(data2);   g_free(data);  
  }
  else
  {
  
      if (tetra)
     {
      gts_isosurface_tetra (surf1, g, fill_function, ddata1, isolevel);
      gts_isosurface_tetra (surf2, g, fill_function, ddata2, isolevel);
      gts_isosurface_tetra (surfg, g, fill_function, ddata,  isolevel);
     }
     else
     {
      gts_isosurface_cartesian (surf1, g, fill_function, ddata1, isolevel);
      gts_isosurface_cartesian (surf2, g, fill_function, ddata2, isolevel);
      gts_isosurface_cartesian (surfg, g, fill_function, ddata,  isolevel);
     }

     /* Free data */
     g_free(ddata1);    g_free(ddata2);   g_free(ddata);    
  }
   
  /* Display summary information about surfaces */
  if (verbose)
  {
    fprintf(stderr,"\n#FLUID 1 surface S1");
    gts_surface_print_stats (surf1, stderr);
    fprintf(stderr,"\n#FLUID 2 surface S2");
    gts_surface_print_stats (surf2, stderr);
    fprintf(stderr,"\n#GRAIN surface SG");
    gts_surface_print_stats (surfg, stderr);
  }
  fflush(stderr);

 /*   No checks (self-intersection etc.) are performed on surfaces. 
  *   They should be consistently oriented and have triangles of good 
  *   quality (since produced by marching cubes algorithm).
  *   When checking for surface triangle intersections, the only thing checked
  *   is that intersecting triangles have the same set of vertices.
  */

  /* Flag vertices bitwise */
  //initialize
  gts_surface_foreach_vertex(surf1,(GtsFunc)flag_initial,&S1);
  gts_surface_foreach_vertex(surf2,(GtsFunc)flag_initial,&S2);
  gts_surface_foreach_vertex(surfg,(GtsFunc)flag_initial,&SG);

  gts_surface_foreach_face(surf1,(GtsFunc)flag_initial,&S1);
  gts_surface_foreach_face(surf2,(GtsFunc)flag_initial,&S2);
  gts_surface_foreach_face(surfg,(GtsFunc)flag_initial,&SG);

  /* Build bounding box trees for surfaces - used for efficient search*/
  tree1 = gts_bb_tree_surface (surf1);
  tree2 = gts_bb_tree_surface (surf2);
  treeg = gts_bb_tree_surface (surfg);

  /* Process triangles with common vertices - this will bitwise
   * flag vertices that belong to more than one surface
  */
  gts_bb_tree_traverse_overlapping (tree1, tree2,
		       (GtsBBTreeTraverseFunc) intersect_triangles_new,NULL);
  gts_bb_tree_traverse_overlapping (treeg, tree1,
		       (GtsBBTreeTraverseFunc) intersect_triangles_new,NULL);
  gts_bb_tree_traverse_overlapping (treeg, tree2,
		       (GtsBBTreeTraverseFunc) intersect_triangles_new,NULL);
  
  /* Now flag triangles according to its vertices */
  gts_surface_foreach_face(surf1,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(surf2,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(surfg,(GtsFunc)flag_triangle_from_vert,NULL);

  /* Get actual surfaces from all the flags*/
  surf12 =  get_surface_from_triangle_flags(surf1,S1+S2,S_DUMMY);  
  surf21 =  get_surface_from_triangle_flags(surf2,S1+S2,S_DUMMY);
  
  surfg1 =  get_surface_from_triangle_flags(surfg,SG+S1,S_DUMMY);
  surf1g =  get_surface_from_triangle_flags(surf1,SG+S1,S_DUMMY);

  surfg2 =  get_surface_from_triangle_flags(surfg,SG+S2,S_DUMMY);
  surf2g =  get_surface_from_triangle_flags(surf2,SG+S2,S_DUMMY);
 
  all = S1 + S2 + SG;
  surf_1only = get_surface_from_triangle_flags(surf1,all,S1);
  surf_2only = get_surface_from_triangle_flags(surf2,all,S2);
  surf_gonly = get_surface_from_triangle_flags(surfg,all,SG);

  
  /* Due to potentially different triangulations in surf1 and surfg
     of the same intersecting patch, e.g.
        ____      and   ____
       |\  |           |  /|
       | \ |           | / |
       |  \|           |/  |
       -----           -----

 surface 1 and g intersections in surf1g and surfg1 might differ.
  */
  //n1g = gts_surface_face_number(surf1g);
  //ng1 = gts_surface_face_number(surfg1);
  //printf("\nsurf1g #face %d  surfg1  #face %d\n",n1g,ng1);
  

  if( color_plot )
  {
     /* Geomview files of (colored) surfaces are output by default.
      * Each file has local axes added so they can be superimposed at will.
      */ 

     //individual surfaces
     gts_surface_foreach_face (surf1, (GtsFunc)set_red_color,NULL);

     sprintf(fname,"%s_fld1.list",base);
     gts_output_surface_oogl_spec(surf1,&g,fname);

     gts_surface_foreach_face (surf2, (GtsFunc)set_green_color, NULL);
     sprintf(fname,"%s_fld2.list",base);
     gts_output_surface_oogl_spec(surf2,&g,fname);

     gts_surface_foreach_face (surfg, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_grain.list",base);
     gts_output_surface_oogl_spec(surfg,&g,fname);

     //intersections file
     sprintf(fname,"%s_intersect.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);

     //gts_surface_write(surf12,stdout);
     Add_surf(surf12,g,fp,"S1 & S2",set_yellow_color,base,"_yellow.list")
     Add_surf(surf1g,g,fp,"S1 & SG",set_magenta_color,base,"_magenta.list")
     Add_surf(surfg2,g,fp,"SG & S2",set_cyan_color,base,"_cyan.list")

     Add_surf(surf_1only,g,fp,"S1 \\ (S2 U SG)",set_red_color,base,"_red.list")
     Add_surf(surf_2only,g,fp,"S2 \\ (S1 U SG)",set_green_color,base,"_green.list")
     Add_surf(surf_gonly,g,fp,"SG \\ (S1 U S2)",set_blue_color,base,"_blue.list")
     fclose (fp);
  }
  
   
  /** Compute areas of all surfaces **/  
  area_12 = gts_surface_area (surf12);
  area_g1 = gts_surface_area (surfg1);
  area_g2 = gts_surface_area (surfg2);

  area_1only = gts_surface_area (surf_1only);
  area_2only = gts_surface_area (surf_2only);
  area_gonly = gts_surface_area (surf_gonly);

   
  sprintf(fname,"%s_areas",base);
  fp = fopen(fname,"a");
  if( web_plot) fprintf(fp,"\nBlob %d",blob_id);

  fprintf(fp,"\nAreas S12 %g S1G %g S2G %g S1only %g S2only %g SGonly %g",
                  area_12, area_g1, area_g2,area_1only,area_2only,area_gonly);
  //Area 1 presumed closed (not checked), so output total area & volume
  area_1 =  area_12 + area_g1 + area_1only;
  vol_1 = gts_surface_volume(surf1);
  fprintf(fp,"\nTotal_S1_area %g  volume %g ratio %g",area_1,
	                                                 vol_1,area_1/vol_1);
  fclose(fp);
  
  /***
  * pointer to angle values for applicable vertices will be stored in 'reserved'
  * so we reset the values beforehand
  * DO NOT USE reserved field for anything else!!!!
  ***/
  gts_surface_foreach_vertex(surf_1only,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(surf_2only,(GtsFunc)gts_object_reset_reserved,NULL);

  
  //initialize Pass_data_vca structure
  vca1.surf_Aonly = surf_1only;
  vca1.surfAB = surf12;
  vca1.surfAC = surf1g;
  vca1.flagAB = S1 + S2;
  vca1.flagAC = S1 + SG;
  vca1.thetaA = NULL; // if used to store a list, has to be initialized to NULL
  
  if( dbg_plot )
  {
     sprintf(fname,"%s_dbg_normals.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);
     vca1.dbg_fp = fp;
  }
  else
  {
    vca1.dbg_fp = NULL;
  }
    
  find_all_vertex_angles(&vca1);
  if( dbg_plot ) fclose(fp);
  
 
  //initialize Pass_data_vca structure
  vca2.surf_Aonly = surf_2only;
  vca2.surfAB = surf21;
  vca2.surfAC = surf2g;
  vca2.flagAB = S2 + S1;
  vca2.flagAC = S2 + SG;
  vca2.thetaA = NULL; 
  vca2.dbg_fp = NULL;
  
  find_all_vertex_angles(&vca2);
   
  master_vert_list = NULL;
  gts_surface_foreach_vertex(surf_1only,(GtsFunc)build_master_list,
                                                       &master_vert_list); 
  gts_surface_foreach_vertex(surf_2only,(GtsFunc)build_master_list,
                                                       &master_vert_list);
   //gather statistics on angles
  sprintf(fname,"%s_cont_angle",base);
  fp  = fopen(fname,"a");
  if( web_plot) fprintf(fp,"\nBlob %d",blob_id);
  get_angle_stats(master_vert_list,fp);
  fclose(fp);
  
  if( print_pointwise )
  {//output all theta1 and theta2 pointwise values to separate files
      FILE *fp1, *fp2;
      
      sprintf(fname,"%s_%d_theta1",base,blob_id);
      fp1 = fopen(fname,"w");
      sprintf(fname,"%s_%d_theta2",base,blob_id);
      fp2 = fopen(fname,"w");
      print_all_angles(master_vert_list,fp1,fp2);
      fclose(fp1);
      fclose(fp2);
  }

  if( web_plot )
  {
     /* Specific plots I need for webpage
      */ 

     //grain surface
     gts_surface_foreach_face (surfg, (GtsFunc)set_no_color, NULL);
     sprintf(fname,"%s_%d_grain.list",base,blob_id);
     gts_output_surface_oogl_spec(surfg,&g,fname);

     //individual blob file with fld-fld and fld-grain surface
     sprintf(fname,"%s_%d.list",base,blob_id);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g.x,g.y,g.z,g.x+g.nx,g.y+g.ny,g.z+g.nz);

     //gts_surface_write(surf12,stdout);
     //this merges is for plot purposes only
     gts_surface_merge(surf12,surf_1only);
     Add_surf_only(surf12,g,fp,"S1 & S2 and S1 \\ (S2 U SG)",set_red_color)
     Add_surf_only(surf1g,g,fp,"S1 & SG",set_no_color)
    
     fclose (fp);
  }
 				
				       
   /* destroy surfaces we don't need*/
  gts_object_destroy (GTS_OBJECT (surf1));
  gts_object_destroy (GTS_OBJECT (surf2));
  gts_object_destroy (GTS_OBJECT (surfg));
  gts_object_destroy (GTS_OBJECT (surf12));
  gts_object_destroy (GTS_OBJECT (surfg1));
  gts_object_destroy (GTS_OBJECT (surf1g));
  gts_object_destroy (GTS_OBJECT (surfg2));
  gts_object_destroy (GTS_OBJECT (surf_1only));
  gts_object_destroy (GTS_OBJECT (surf_2only));
  gts_object_destroy (GTS_OBJECT (surf_gonly));
  
    /*Destroy trees & bounding boxes */
  gts_bb_tree_destroy (tree1, TRUE);
  gts_bb_tree_destroy (tree2, TRUE);
  gts_bb_tree_destroy (treeg, TRUE);
  
  MEM_STATS();
							 
  return 0;
}


/* void fill_from_data()
*  Fills 2D double array 'f' from packed unsigned char array 'data'

*/
static void fill_from_data(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data)
{
  guint   i, j, nxy = g.nx * g.ny;
  guchar  *pdata;
  
  pdata  = data + k*nxy;
  for (j = 0; j < g.ny; j++)
      for (i = 0; i < g.nx; i++, pdata++)
               f[i][j] = *pdata ;
}

/* 
*  void fill_from_double_data()
*  Fills 2D double array 'f' from packed unsigned char array 'data'
*  Definition of this function recquired in order to use
*  marching cubes isosurfacing functions
*  gts_isosurface_cartesian() and gts_isosurface_tetra()
*/

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




void  create_test_example(gdouble **pdata1, gdouble **pdata2, gdouble **pdata,Cube_Info *cube)
{
    gint     n = 2;                        //building block size
    gint     nx = 5,     ny = 8,   nz = 2; //number of building blocks
    gint     cutx = 3, cuty = 4;
    
    gdouble  *data1, *data2, *data;
    
    guint    i, j, k, ind;

    cube->xs = cube->ys = cube->zs = 0;
    cube->nx = n * nx;
    cube->ny = n * ny;
    cube->nz = n * nz;

    cube->nxy  = cube->nx  * cube->ny;
    cube->nxyz = cube->nxy * cube->nz;

    data1 = (gdouble *)g_malloc(cube->nxyz * sizeof(gdouble));
    data2 = (gdouble *)g_malloc(cube->nxyz * sizeof(gdouble));
    data =  (gdouble *)g_malloc(cube->nxyz * sizeof(gdouble));
    
    ind = 0;
    for(k = 0; k < cube->nz; k++)
       for(j = 0; j < cube->ny; j++)
          for(i = 0; i < cube->nx; i++,ind++)
	    {
	        if( (i >= cutx*n) && (j < cuty*n) ) data1[ind] = 0;
                else                                data1[ind] = 1;

                if( (i >= cutx*n) && (j >= cuty*n) ) data2[ind] = 0;
                else                                 data2[ind] = 1;
		
		if( data1[ind] == 1 && data2[ind] == 1 ) data[ind] = 1;
		else                                     data[ind] = 0;
            }

    *pdata1 = data1;
    *pdata2 = data2;
    *pdata = data;
}


