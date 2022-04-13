/* GTS - Library for the manipulation of triangulated surfaces
 * Copyright (C) 1999 Stéphane Popinet
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "gts.h"
#include "iso_modified.h"

static guint32 all, S1 = 1, S2 = 2, SG = 4, S0 = 0; //bitwise surface flags

/* coordinates of the edges of the cube (see doc/isocube.fig) */
static guint c[12][4] = {
  {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 1, 0},
  {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 1, 0, 1}, {1, 1, 0, 0},
  {2, 0, 0, 0}, {2, 1, 0, 0}, {2, 1, 1, 0}, {2, 0, 1, 0}};

/* first index is the edge number, second index is the edge orientation 
   (RIGHT or LEFT), third index are the edges which this edge may connect to
   in order */
static guint edge[12][2][3] = {
  {{9, 1, 8}, {4, 3, 7}},   /* 0 */
  {{6, 2, 5}, {8, 0, 9}},   /* 1 */
  {{10, 3, 11}, {5, 1, 6}}, /* 2 */
  {{7, 0, 4}, {11, 2, 10}}, /* 3 */
  {{3, 7, 0}, {8, 5, 11}},  /* 4 */
  {{11, 4, 8}, {1, 6, 2}},  /* 5 */
  {{2, 5, 1}, {9, 7, 10}},  /* 6 */
  {{10, 6, 9}, {0, 4, 3}},  /* 7 */
  {{5, 11, 4}, {0, 9, 1}},  /* 8 */
  {{1, 8, 0}, {7, 10, 6}},  /* 9 */
  {{6, 9, 7}, {3, 11, 2}},  /* 10 */
  {{2, 10, 3}, {4, 8, 5}}   /* 11 */
};

static void ** malloc2D (guint nx, guint ny, gulong size)
{
  void ** m = g_malloc (nx*sizeof (void *));
  guint i;

  for (i = 0; i < nx; i++)
    m[i] = g_malloc0 (ny*size);

  return m;
}

static void free2D (void ** m, guint nx)
{
  guint i;

  g_return_if_fail (m != NULL);

  for (i = 0; i < nx; i++)
    g_free (m[i]);
  g_free (m);
}



/**
 * gts_iso_slice_modified_new:
 * @nx: number of vertices in the x direction.
 * @ny: number of vertices in the y direction.
 *
 * Returns: a new #GtsIsoSliceModified.
 */
GtsIsoSliceModified * gts_iso_slice_modified_new (guint nx, guint ny)
{
  GtsIsoSliceModified * slice;

  g_return_val_if_fail (nx > 1, NULL);
  g_return_val_if_fail (ny > 1, NULL);

  slice = g_malloc (sizeof (GtsIsoSliceModified));

  slice->vertices = g_malloc (3*sizeof (OrientedVertexModified **));
  slice->vertices[0] = 
    (OrientedVertexModified **) malloc2D (nx, ny, sizeof (OrientedVertexModified));
  slice->vertices[1] = 
    (OrientedVertexModified **) malloc2D (nx - 1, ny, sizeof (OrientedVertexModified));
  slice->vertices[2] = 
    (OrientedVertexModified **) malloc2D (nx, ny - 1, sizeof (OrientedVertexModified));
  slice->nx = nx;
  slice->ny = ny;

  return slice;
}


#define ORDER1 0
#define ORDER0 1
 
/**
 * gts_iso_slice_fill_cartesian_modified:
 * @slice: a #GtsIsoSliceModified.
 * @g: a #GtsCartesianGrid.
 * @f1: values of the function for plane z = @g.z.
 * @f2: values of the function for plane z = @g.z + @g.dz.
 * @iso: isosurface value.
 * @klass: a #GtsVertexClass.
 *
 * Fill @slice with the coordinates of the vertices defined by 
 * f1 (x,y,z) = @iso and f2 (x, y, z) = @iso.
 */
 
 /* f < iso is assumed to be where the fluid 
   we (linearly) extrapolate fluid values (stored in ff* )
   from f < iso phase onto the f=iso surface vertices
   
 */
void gts_iso_slice_fill_cartesian_modified (GtsIsoSliceModified * slice,
				   GtsCartesianGrid g,
				   gdouble ** f_1,
				   gdouble ** f0,
				   gdouble ** f1,
				   gdouble ** f2,
				   gdouble ** f3,
				   gdouble ** f4,
				   gdouble ** ff_1,
				   gdouble ** ff0,
				   gdouble ** ff1,
				   gdouble ** ff2,
				   gdouble ** ff3,
				   gdouble ** ff4,
				   gdouble iso,
				   gdouble iso_fluid,
				   GtsVertexClass * klassG)
{
  OrientedVertexModified *** vertices;
  guint i, j;
  gdouble x, y, rho, one_minus_rho;
  
  gdouble a,b,c, vplus1, vplus2, fluid_value, test;
  guint32 flag;

  g_return_if_fail (slice != NULL);
  g_return_if_fail (f1 != NULL);
  g_return_if_fail (ff1 != NULL);
  
  vertices = slice->vertices;

  //printf("\nz %g\n",g.z);
#define Z  
#ifdef Z  
  if (f2 && ff2)
    for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
      {
	gdouble v1 = f1[i][j] - iso;
	gdouble v2 = f2[i][j] - iso;
	if ((v1 >= 0. && v2 < 0.) || (v1 < 0. && v2 >= 0.))
	{
	  rho = v1/(v1 - v2);
	  vertices[0][i][j].v = 
	    gts_vertex_new (klassG,
			    x, y, g.z + g.dz*rho);
	    vertices[0][i][j].orientation = v2 >= 0. ? RIGHT : LEFT;
	  
	  
	  //printf("\n i %d j %d, v1 %g v2 %g rho %g",i,j,v1,v2,rho);
	  
	  /* figure out fluid value for the vertex */
	  if ( vertices[0][i][j].orientation == RIGHT ) 
	      /* v1 is the fluid vertex (v2 is grain )
	         find  fluid value of the surface vertex
		 by extrapolation from slices below */
	  {
	     a = ff1[i][j]; /*this is the simplest extrapol. value we could use */
	     //printf("  Z ff1 %g",a);
	     /* check out if there are more fluid vertices to use in extrapolation */
	     if (f0 && ff0 )  vplus1 = f0[i][j] - iso;
	     else             vplus1 = 0;
	     
	     if (f_1 && ff_1) vplus2 = f_1[i][j] - iso;
	     else             vplus2 = 0;
	     
	     if (ORDER0) vplus1 = 0; 
	     if (ORDER1 || ORDER0) vplus2 = 0;	     	     	     
	     if ( vplus1 < 0. && vplus2 < 0.) 
	     {	     	     	        
	      
		   b = -1.5*a + 2.0*ff0[i][j] - 0.5*ff_1[i][j];
		   c =  0.5*a  -    ff0[i][j] + 0.5*ff_1[i][j];
		   fluid_value = a + b * rho + c*rho*rho;
             }
	     else if ( vplus1 < 0.  )
             {
		   b = ff0[i][j] - a;
		   fluid_value = a + b * rho; 
	     }
	     else
	     {
	        fluid_value = a;
	     }	     	     	     	     	  	  
	  }
	  else  /* v2 is the fluid vertex
	        find  fluid value of the surface vertex
		by extrapolation from slices below */
	  {
	     a = ff2[i][j];
	     //printf("  Z ff2 %g",a);
	     one_minus_rho = 1 - rho;
	     
	     if (f3 && ff3)  vplus1 = f3[i][j] - iso;
	     else            vplus1 = 0;
	     	     
	     if (f4 && ff4)  vplus2 = f4[i][j] - iso;
	     else            vplus2 = 0;
	     
	     if (ORDER0)  vplus1 = 0; 
	     if (ORDER1 || ORDER0) vplus2 = 0;	     	     	     
	     if ( vplus1 < 0. && vplus2 < 0. ) 
	     {
		   b = -1.5*a + 2.0*ff3[i][j] - 0.5*ff4[i][j];
		   c =  0.5*a -     ff3[i][j] + 0.5*ff4[i][j];
		   
		   fluid_value = a + b * one_minus_rho + c*one_minus_rho*one_minus_rho;
             }
	     else if ( vplus1 < 0. )
	     {
		   b = ff3[i][j] - a;
		   fluid_value = a + b * one_minus_rho;
	     }
	     else
	     {
	        fluid_value = a;
	     }	        	     	     	  	  
	  }
	  
	  vertices[0][i][j].fluid_value = fluid_value;
	  
	  test = fluid_value - iso_fluid;
	  if( test < 0. )
	    flag = S1 + SG;
	  else if ( test > 0. )
	    flag = S2 + SG;
	  else
	    flag = SG;  
	    
	  GTS_OBJECT(vertices[0][i][j].v)->flags = flag;
	  //printf("  flag %d",flag);
	}
	else
	  vertices[0][i][j].v = NULL;
      }
 #endif     

#define X
#ifdef X      
  for (i = 0, x = g.x; i < g.nx - 1; i++, x += g.dx)
    for (j = 0, y = g.y; j < g.ny; j++, y += g.dy) 
    {
      gdouble v1 = f1[i][j] - iso;
      gdouble v2 = f1[i+1][j] - iso;
      if ((v1 >= 0. && v2 < 0.) || (v1 < 0. && v2 >= 0.))
      {
        rho = v1/(v1 - v2);
	vertices[1][i][j].v = 
	  gts_vertex_new (klassG, x + g.dx*rho, y, g.z);
	vertices[1][i][j].orientation = v2 >= 0. ? RIGHT : LEFT;
	
        //printf("\n i %d j %d, v1 %g v2 %g rho %g",i,j,v1,v2,rho);
        
	
	/* figure out fluid value for the vertex */	
	 if ( vertices[1][i][j].orientation == RIGHT )  
	        /* v1 is the fluid vertex (v2 is grain )
	         find  fluid value of the surface vertex
		 by extrapolation "lesser" x-values */
	  {
	     a = ff1[i][j];
	     //printf("  X ff1 %g",a);

	     if (i >= 1)  vplus1 = f1[i-1][j] - iso;
	     else        vplus1 = 0;
	     	
	     if (i >= 2)  vplus2 = f1[i-2][j] - iso;
	     else        vplus2 = 0;
		
	     if (ORDER0)  vplus1 = 0; 
	     if (ORDER1 || ORDER0) vplus2 = 0;		     	     	     
	     if ( vplus1 < 0. && vplus2 < 0. ) 
	     {
		   b = -1.5*a + 2.0*ff1[i-1][j] - 0.5*ff1[i-2][j];
		   c =  0.5*a -     ff1[i-1][j] + 0.5*ff1[i-2][j]; 
		   fluid_value = a + b * rho + c*rho*rho;
	     }
             else if ( vplus1 < 0.)
	     {
		   b = ff1[i-1][j] - a;
		   fluid_value = a + b * rho;
	     } 
	     else
	     {
	        fluid_value = a;
	     }
	         	     	     	  	  
	  }
	  else  /* v2 is the fluid vertex,
	         find  fluid value of the surface vertex
		 by extrapolation  "upper" x-values */
	  {
	     a = ff1[i+1][j];
	     one_minus_rho = 1 - rho;	     
	    
	     //printf("  X ff1(i+1) %g",a);
	     
	     if (i < g.nx-2)  vplus1 = f1[i+2][j] - iso;
	     else             vplus1 = 0;
	     
	     if (i < g.nx-3)  vplus2 = f1[i+3][j] - iso;
	     else             vplus2 = 0;
	     	 
	     if (ORDER0)  vplus1 = 0;
	     if (ORDER1 || ORDER0) vplus2 = 0;	     	     	     
	     if ( vplus1 < 0. && vplus2 < 0. ) 
	     {
		   b = -1.5*a + 2.0*ff1[i+2][j] - 0.5*ff1[i+3][j];
		   c =  0.5*a -     ff1[i+2][j] + 0.5*ff1[i+3][j];
		   fluid_value = a + b * one_minus_rho + c*one_minus_rho*one_minus_rho;
	     }
	     else if ( vplus1 < 0. )
	     {
		   b = ff1[i+2][j] - a;
		   fluid_value = a + b * one_minus_rho;
	     } 
	     else
	     {
	        fluid_value = a;
	     }	    	     	     	     	  	  
	  }	 


	  vertices[1][i][j].fluid_value = fluid_value;
	  
	  test = fluid_value - iso_fluid;
	  if( test < 0. )
	    flag = S1 + SG;
	  else if ( test > 0. )
	    flag = S2 + SG;
	  else
	    flag = SG;  
	    
	  GTS_OBJECT(vertices[1][i][j].v)->flags = flag;
	  //printf("  flag %d",flag);	
      }
      else
	vertices[1][i][j].v = NULL;
    }
#endif    
    

#define Y
#ifdef  Y    
  for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
    for (j = 0, y = g.y; j < g.ny - 1; j++, y += g.dy)
    {
      gdouble v1 = f1[i][j] - iso;
      gdouble v2 = f1[i][j+1] - iso;
      if ((v1 >= 0. && v2 < 0.) || (v1 < 0. && v2 >= 0.))
      {
        rho = v1/(v1-v2);
	vertices[2][i][j].v = 
	  gts_vertex_new (klassG, x, y + g.dy*rho, g.z);
	vertices[2][i][j].orientation = v2 >= 0. ? RIGHT : LEFT;
	
        
	//printf("\n i %d j %d, v1 %g v2 %g rho %g",i,j,v1,v2,rho);
		
	/* figure out fluid value for the vertex */	
	 if ( vertices[2][i][j].orientation == RIGHT ) 
	            /* v1 is the fluid vertex
		      extrapolation from "lesser" y-values */
	  {
	     a = ff1[i][j];	     
	     //printf("  Y ff1 %g",a);
	     if (j >= 1)  vplus1 = f1[i][j-1] - iso;
	     else        vplus1 = 0;
	     
	     if (j >= 2)  vplus2 = f1[i][j-2] - iso;
	     else        vplus2 = 0;
	     
	     if (ORDER0)  vplus1 = 0;
	    if (ORDER1 || ORDER0) vplus2 = 0;	     	     	     
	     if ( vplus1 < 0. && vplus2 < 0. )
	     {
		   b = -1.5*a + 2.0*ff1[i][j-1] - 0.5*ff1[i][j-2];
		   c =  0.5*a -     ff1[i][j-1] + 0.5*ff1[i][j-2]; 
		   fluid_value = a + b * rho + c*rho*rho;
             }
	     else if ( vplus1 < 0. )
	     {
		   b = ff1[i][j-1] - a;
		   fluid_value = a + b * rho;
             } 	     
	     else
	     {
	        fluid_value = a;
	     }
	        	     	     	  	  
	  }
	  else  /* v2 is the fluid vertex
	         extrapolation from "upper" y-values */
	  {
	     a = ff1[i][j+1];
	     one_minus_rho = 1 -rho;
	     
	     //printf("  Y ff1(j+1) %g",a);
	     if (j < g.ny-2)  vplus1 = f1[i][j+2] - iso;
	     else             vplus1 = 0;
	     
	     if (j < g.ny-3)  vplus2 = f1[i][j+3] - iso;
	     else             vplus2 = 0;
	     
	     	     
	     if (ORDER0)  vplus1 = 0;
	     if (ORDER1 || ORDER0) vplus2 = 0;	 	     	     	     
	     if ( vplus1 < 0. && vplus2 < 0. )
	     {
		   b = -1.5*a + 2.0*ff1[i][j+2] - 0.5*ff1[i][j+3];
		   c =  0.5*a -     ff1[i][j+2] + 0.5*ff1[i][j+3];
		   fluid_value = a + b * one_minus_rho + c*one_minus_rho*one_minus_rho;
	     }
	     else if ( vplus1 < 0. )
	     {
		   b = ff1[i][j+2] - a;
		   fluid_value =  a + b * one_minus_rho;
	     }
	     else
	     {
	           fluid_value = a;
	     }	     	     	     	     	  	  
	  }
	 	
	  vertices[2][i][j].fluid_value = fluid_value;
	  
	  test = fluid_value - iso_fluid;
	   if( test < 0. )
	    flag = S1 + SG;
	  else if ( test > 0. )
	    flag = S2 + SG;
	  else
	    flag = SG;  
	    
	  GTS_OBJECT(vertices[2][i][j].v)->flags = flag;
	  //printf("  flag %d",flag);	  
      }
      else
	vertices[2][i][j].v = NULL;
    }
#endif    
    
}

/**
 * gts_iso_slice_modified_destroy:
 * @slice: a #GtsIsoSliceModified.
 *
 * Free all memory allocated for @slice.
 */
void gts_iso_slice_modified_destroy (GtsIsoSliceModified * slice)
{
  g_return_if_fail (slice != NULL);

  free2D ((void **) slice->vertices[0], slice->nx);
  free2D ((void **) slice->vertices[1], slice->nx - 1);
  free2D ((void **) slice->vertices[2], slice->nx);  
  g_free (slice->vertices);
  g_free (slice);
}


/**
 * gts_isosurface_slice_modified:
 * @slice1: a #GtsIsoSliceModified.
 * @slice2: another #GtsIsoSliceModified.
 * @surface: a #GtsSurface.
 *
 * Given two successive slices @slice1 and @slice2 link their vertices with
 * segments and triangles which are added to @surface.
 */
void gts_isosurface_slice_modified (GtsIsoSliceModified * slice1,
			            GtsIsoSliceModified * slice2,
			            GtsSurface * surface)
{
  guint j, k, l, nx, ny;
  OrientedVertexModified *** vertices[2];
  GtsVertex * va[12];

  g_return_if_fail (slice1 != NULL);
  g_return_if_fail (slice2 != NULL);
  g_return_if_fail (surface != NULL);
  g_return_if_fail (slice1->nx == slice2->nx && slice1->ny == slice2->ny);

  vertices[0] = slice1->vertices;
  vertices[1] = slice2->vertices;
  nx = slice1->nx;
  ny = slice1->ny;

  /* link vertices with segments and triangles */
  for (j = 0; j < nx - 1; j++)
    for (k = 0; k < ny - 1; k++) 
    {
      gboolean cube_is_cut = FALSE;
      
      for (l = 0; l < 12; l++) /* traverse cube edges */
      {
	guint nv = 0, e = l;
	OrientedVertexModified ov = 
	  vertices[c[e][1]][c[e][0]][j + c[e][2]][k + c[e][3]];
	while (ov.v && !GTS_OBJECT (ov.v)->reserved)
	{
	  guint m = 0, * ne = edge[e][ov.orientation];
	  va[nv] = ov.v;
	  nv++;
	  GTS_OBJECT (ov.v)->reserved = surface;
	  ov.v = NULL;
	  while (m < 3 && !ov.v) 
	  {
	    e = ne[m++];
	    ov = vertices[c[e][1]][c[e][0]][j + c[e][2]][k + c[e][3]];
	  }
	}
	/* create edges and faces */
	if (nv > 2) 
	{
	  GtsEdge * e1, * e2, * e3;
	  guint m;
	  if (!(e1 = GTS_EDGE (gts_vertices_are_connected (va[0], va[1]))))
	  {
	    e1 = gts_edge_new (surface->edge_class, va[0], va[1]);
	  }  
	  for (m = 1; m < nv - 1; m++) 
	  {
	    if (!(e2 = GTS_EDGE (gts_vertices_are_connected (va[m], va[m+1]))))
	    {
	      e2 = gts_edge_new (surface->edge_class, va[m], va[m+1]);
	    }
	    if (!(e3 = GTS_EDGE (gts_vertices_are_connected (va[m+1], va[0]))))
	    {
	      e3 = gts_edge_new (surface->edge_class, va[m+1], va[0]);
	    }
	    gts_surface_add_face (surface, 
				  gts_face_new (surface->face_class,
						e1, e2, e3));
	    e1 = e3;
	  }
	}
	if (nv > 0)
	  cube_is_cut = TRUE;
      }
      if (cube_is_cut)
	for (l = 0; l < 12; l++) {
	  GtsVertex * v = 
	    vertices[c[l][1]][c[l][0]][j + c[l][2]][k + c[l][3]].v;
	  if (v)
	    GTS_OBJECT (v)->reserved = NULL;
	}
    }
}

#define SWAP(s1, s2, tmp) (tmp = s1, s1 = s2, s2 = tmp)

#define SWAP_MODIFIED(s_1,s0,s1,s2,s3,s4,tmp) \
    (tmp = s_1, s_1 = s0,s0 = s1, s1 = s2, s2 = s3, s3 = s4, s4 = tmp)

/**
 * gts_isosurface_cartesian_modified:
 * @surface: a #GtsSurface.
 * @g: a #GtsCartesianGrid.
 * @f: a #GtsIsoCartesianFunc.
 * @data: user data to be passed to @f.
 * @iso: isosurface value.
 *
 * Adds to @surface new faces defining the isosurface f(x,y,z) = @iso. By
 * convention, the normals to the surface are pointing toward the positive
 * values of f(x,y,z) - @iso.
 * 
 * In addition, we assume that f(x,y,z) - @iso is "fluid" space that has
 * two fluid phases. Fluid phase data_fluid is well defined only within
 * fluid phase and vertices of the isosurface will have fluid data
 * extrapolated from those values.
 *
 * The user function @f is called successively for each value of the z 
 * coordinate defined by @g. It must fill the corresponding (x,y) plane with
 * the values of the function for which the isosurface is to be computed.
 */
void gts_isosurface_cartesian_modified (GtsSurface * surfaceG,
			       GtsCartesianGrid g,
			       GtsIsoCartesianFunc f,
			       gpointer data,
			       gpointer data_fluid,
			       gdouble iso,
			       gdouble iso_fluid)
{
  gdouble **tmp;
  gdouble **f_1,  **f0,  **f1,  **f2,  **f3,  **f4;
  gdouble **ff_1, **ff0, **ff1, **ff2, **ff3, **ff4;
  
  gdouble **pf_1,  **pf0,  **pf1,  **pf2,  **pf3,  **pf4;
  gdouble **pff_1, **pff0, **pff1, **pff2, **pff3, **pff4;
  
  GtsIsoSliceModified * slice1, * slice2, *tmp_slice;
  
  guint i;

  g_return_if_fail (surfaceG != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (g.nx > 1);
  g_return_if_fail (g.ny > 1);
  g_return_if_fail (g.nz > 1);

  /* allocate space for two slices */
  slice1 = gts_iso_slice_modified_new (g.nx, g.ny);
  slice2 = gts_iso_slice_modified_new (g.nx, g.ny);
  
  f_1 = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  f0  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  f1  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  f2  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  f3  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  f4  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
   
  ff_1 = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  ff0  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  ff1  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  ff2  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  ff3  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  ff4  = (gdouble **) malloc2D (g.nx, g.ny, sizeof (gdouble));
  

  pf_1 = f_1; pf0 = f0; pf1 = f1; pf2 = f2; pf3 = f3; pf4 = f4;
  pff_1 = ff_1; pff0 = ff0; pff1 = ff1; pff2 = ff2; pff3 = ff3; pff4 = ff4;
  
  (*f) (pf_1, g, 0, data);
  (*f) (pff_1, g, 0, data_fluid);
  
  (*f) (pf0, g, 0, data);
  (*f) (pff0, g, 0, data_fluid);
  
  (*f) (pf1, g, 0, data);
  (*f) (pff1, g, 0, data_fluid);

  g.z += g.dz;  
  (*f) (pf2, g, 1, data);
  (*f) (pff2, g, 1, data_fluid);
  
  if( g.nz >= 3 )
  {       
     g.z += g.dz;  
     (*f) (pf3, g, 2, data);
     (*f) (pff3, g, 2, data_fluid);
     
     if( g.nz >= 4 )
     { 
       g.z += g.dz;  
       (*f) (pf4, g, 3, data);
       (*f) (pff4, g, 3, data_fluid);
       
       g.z -=g.dz;
     }
     else
     {
       /* fill data from from position 2 */
       (*f) (pf4, g, 2, data);
       (*f) (pff4, g, 2, data_fluid);
     }     
     g.z -=g.dz;
  }
  else
  {  /* fill data from from  position 1 */
     (*f) (pf3, g, 1, data);
     (*f) (pff3, g, 1, data_fluid);
     
     (*f) (pf4, g, 1, data);
     (*f) (pff4, g, 1, data_fluid);
  }
   
  g.z -= g.dz;
  
  gts_iso_slice_fill_cartesian_modified (slice1, g, 
                                pf_1,pf0,pf1, pf2,pf3,pf4,
				pff_1,pff0,pff1,pff2,pff3,pff4,
				iso,iso_fluid, 
				surfaceG->vertex_class);
  g.z += g.dz;
  
  for (i = 2; i < g.nz; i++) 
  {
    if ( i+2 < g.nz )
    {  
       g.z += 2*g.dz;    
    
      (*f) (pf_1, g, i+2, data); 
      tmp = pf_1; pf_1 = pf0; pf0 = pf1; pf1 = pf2; pf2 = pf3; pf3 = pf4; pf4 = tmp;
    			      
      (*f) (pff_1, g, i+2, data_fluid);
      tmp = pff_1; pff_1 = pff0; pff0 = pff1;pff1 = pff2; pff2 = pff3; pff3 = ff4;
      ff4 = tmp;
      g.z -= 2*g.dz;
    }
   else  if ( i+1 < g.nz )
    {
        pf_1 = pf0; pf0 = pf1; pf1 = pf2; pf2 = pf3; pf3 = pf4;
	pff_1 = pff0; pff0 = pff1; pff1 = pff2; pff2 = pff3; pff3 = pff4;
    }
    else 
    {
        pf_1 = pf0; pf0 = pf1; pf1 = pf2; pf2 = pf3;
	pff_1 = pff0; pff0 = pff1; pff1 = pff2; pff2 = pff3;
    }
    
    
    gts_iso_slice_fill_cartesian_modified (slice2,  g, 
                                pf_1,pf0,pf1, pf2,pf3,pf4,
				pff_1,pff0,pff1,pff2,pff3,pff4,
				iso,iso_fluid, 
				surfaceG->vertex_class);
    g.z += g.dz;    
    gts_isosurface_slice_modified (slice1, slice2, surfaceG);
    SWAP (slice1, slice2, tmp_slice);
  }  
  				

  /* free memory */
  
  gts_iso_slice_modified_destroy (slice1);
  gts_iso_slice_modified_destroy (slice2);
 
 /* swapping above does some strange business - segfault  
  free2D ((void **) f_1, g.nx);  
  free2D ((void **) f0, g.nx);  
  free2D ((void **) f1, g.nx);  
  free2D ((void **) f2, g.nx);
  free2D ((void **) f3, g.nx);  
  free2D ((void **) f4, g.nx);
  
  
  free2D ((void **) ff_1, g.nx);
  free2D ((void **) ff0, g.nx);
  free2D ((void **) ff1, g.nx);
  free2D ((void **) ff2, g.nx);
  free2D ((void **) ff3, g.nx);
  free2D ((void **) ff4, g.nx); */
}
