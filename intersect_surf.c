#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <float.h>

#include "config.h"
#include "gts.h"
#include "intersect_surf.h"
#include "mem.h"
#include "viz_util.h"

/* 
*  void fill_from_data()
*  Fills 2D double array 'f' from packed unsigned char array 'data'
*  Definition of this function recquired in order to use
*  marching cubes isosurfacing functions
*  gts_isosurface_cartesian() and gts_isosurface_tetra()
*/

static void fill_from_data(
      gdouble          **f,
      GtsCartesianGrid   g,
      guint              k,
      gpointer           data)
{
  guint	i, j;
  guchar  *pdata;
  gint     nxy = g.nx * g.ny;
  
  pdata = data + k*nxy;
  for( j = 0; j < g.ny; j++ )
	for( i = 0; i < g.nx; i++, pdata++ )
		f[i][j] = (gdouble)*pdata ;
}


/*
* void   intersect_surfaces_from_data_base()
* Returns pointer (SurfacePack *)pack_ptr
*   pack_ptr->surf[1,2,g] is surface of data[1,2,g], and all other surfaces
*   in *packptr are their mutual intersections
*/

SurfacePack *intersect_surfaces_from_data_base(
       guchar     *data1, //assumed 0 where fluid1 is, otherwise 1
       guchar     *data2, //assumed 0 where fluid2 is, otherwise 1
       guchar     *data,  //assumed 0 where grain  is, otherwise 1
       GtsCartesianGrid *g, 
       guint32    S1,  //bitwise mark surf1
       guint32    S2,  //bitwise mark surf2
       guint32    SG,  //bitwise mark surfg
       gboolean   tetra,
       gboolean   verbose)	
{
  gdouble isolevel;
  GtsIsoCartesianFunc func = fill_from_data;

  SurfacePack  *pack;
  GNode	*tree1, *tree2, *treeg;

  guint32  all, S_DUMMY = 0;
  
  //allocate memory
  pack = (SurfacePack *)MALLOC(SPsz);

 /* Initialize surface structures */
  pack->surf1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());
  pack->surf2 = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());
  pack->surfg = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());

  /*
   * Produce triangulated surfaces by marching cubes algorithm.
   * Note that each triangle is oriented so that its normal by default points
   * to larger data values (i.e. 1 in our case).
   * Therefore pack->surf1 normals point outwards from fluid 1 (pack->surf2, pack->surfg similarly).
   */
  isolevel = 0.5;
  
  if( tetra )
  {
    gts_isosurface_tetra (pack->surf1, *g, func, data1, isolevel);
    gts_isosurface_tetra (pack->surf2, *g, func, data2, isolevel);
    gts_isosurface_tetra (pack->surfg, *g, func, data,  isolevel);
  }
  else
  {
    gts_isosurface_cartesian (pack->surf1, *g, func, data1, isolevel);
    gts_isosurface_cartesian (pack->surf2, *g, func, data2, isolevel);
    gts_isosurface_cartesian (pack->surfg, *g, func, data,  isolevel);
  }

  /* Display summary information about surfaces */
  if( verbose )
  {
    fprintf(stderr,"\n#FLUID 1 surface S1");
    gts_surface_print_stats (pack->surf1, stderr);
    fprintf(stderr,"\n#FLUID 2 surface S2");
    gts_surface_print_stats (pack->surf2, stderr);
    fprintf(stderr,"\n#GRAIN surface SG");
    gts_surface_print_stats (pack->surfg, stderr);
  }


 /*   No checks (self-intersection etc.) are performed on surfaces.
  *   They should be consistently oriented and have triangles of good
  *   quality (since produced by marching cubes algorithm).
  *   When checking for surface triangle intersections, the only thing checked
  *   is that intersecting triangles have the same set of vertices.
  */

  /* Flag vertices bitwise */
  //initialize
  gts_surface_foreach_vertex(pack->surf1,(GtsFunc)flag_initial,&S1);
  gts_surface_foreach_vertex(pack->surf2,(GtsFunc)flag_initial,&S2);
  gts_surface_foreach_vertex(pack->surfg,(GtsFunc)flag_initial,&SG);

  gts_surface_foreach_face(pack->surf1,(GtsFunc)flag_initial,&S1);
  gts_surface_foreach_face(pack->surf2,(GtsFunc)flag_initial,&S2);
  gts_surface_foreach_face(pack->surfg,(GtsFunc)flag_initial,&SG);

  /* Build bounding box trees for surfaces - used for efficient search*/
  tree1 = gts_bb_tree_surface (pack->surf1);
  tree2 = gts_bb_tree_surface (pack->surf2);
  treeg = gts_bb_tree_surface (pack->surfg);

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
  gts_surface_foreach_face(pack->surf1,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(pack->surf2,(GtsFunc)flag_triangle_from_vert,NULL);
  gts_surface_foreach_face(pack->surfg,(GtsFunc)flag_triangle_from_vert,NULL);

  /* Get actual surfaces from all the flags*/
  pack->surf12 =  get_surface_from_triangle_flags(pack->surf1,S1+S2,S_DUMMY);
  pack->surf21 =  get_surface_from_triangle_flags(pack->surf2,S1+S2,S_DUMMY);

  pack->surfg1 =  get_surface_from_triangle_flags(pack->surfg,SG+S1,S_DUMMY);
  pack->surf1g =  get_surface_from_triangle_flags(pack->surf1,SG+S1,S_DUMMY);

  pack->surfg2 =  get_surface_from_triangle_flags(pack->surfg,SG+S2,S_DUMMY);
  pack->surf2g =  get_surface_from_triangle_flags(pack->surf2,SG+S2,S_DUMMY);

  all = S1 + S2 + SG;
  pack->surf_1only = get_surface_from_triangle_flags(pack->surf1,all,S_DUMMY);
  pack->surf_2only = get_surface_from_triangle_flags(pack->surf2,all,S_DUMMY);
  pack->surf_gonly = get_surface_from_triangle_flags(pack->surfg,all,S_DUMMY);
  
  /* Display summary information about intersecting surfaces */
  if( verbose )
  {
     /*Note: 
    surf1g is a descendant of surf1, and surfg1 is descendant from surfg
    though they are both intersections of surf and surfg.

    Due to potentially different triangulations in surf1 and surfg
      of the same intersecting patch, e.g.
         ____      and   ____
	|\  |           |  /|
	| \ |           | / |
	|  \|           |/  |
	-----           -----

    surface 1 and g intersections in surf1g and surfg1 might differ.
   */
  
    fprintf(stderr,"\n#FLUID1 and FLUID2 intersection - S1 & S2");
    gts_surface_print_stats (pack->surf12, stderr);
    fprintf(stderr,"\n#FLUID2 and FLUID1 intersection - S2 & S1");
    gts_surface_print_stats (pack->surf21, stderr);
    
    fprintf(stderr,"\n#FLUID1 and GRAIN intersection - S1 & SG");
    gts_surface_print_stats (pack->surf1g, stderr);
    fprintf(stderr,"\n#GRAIN and FLUID1 intersection - SG & S1");
    gts_surface_print_stats (pack->surfg1, stderr);
    
    fprintf(stderr,"\n#FLUID2 and GRAIN intersection - S2 & SG");
    gts_surface_print_stats (pack->surf2g, stderr);
    fprintf(stderr,"\n#GRAIN and FLUID2 intersection - SG & S2");
    gts_surface_print_stats (pack->surfg2, stderr);
    
    fprintf(stderr,"\nFLUID1 only - not intersecting S2 or SG");
    gts_surface_print_stats (pack->surf_1only, stderr);
    fprintf(stderr,"\nFLUID2 only - not intersecting S1 or SG");
    gts_surface_print_stats (pack->surf_2only, stderr);
    fprintf(stderr,"\nGRAIN  only - not intersecting S1 or S2");
    gts_surface_print_stats (pack->surf_gonly, stderr);
  }

   /*Destroy trees & bounding boxes */
  gts_bb_tree_destroy (tree1, TRUE);
  gts_bb_tree_destroy (tree2, TRUE);
  gts_bb_tree_destroy (treeg, TRUE);

  return pack;
}



void	intersect_triangles_new(GtsBBox * bb1, GtsBBox * bb2)
{
  GtsTriangle * t1 = GTS_TRIANGLE (bb1->bounded);
  GtsTriangle * t2 = GTS_TRIANGLE (bb2->bounded);

  match_triangle_vertices(t1,t2);
}



inline gint	match_triangle_vertices(
		GtsTriangle	*t1,
		GtsTriangle	*t2)
{
	gint	equal[3];
	gint	i , j;
	guint32 Fi, Fj;

	GtsVertex *v[3], *w[3];
	GtsPoint  *p1, *p2;
	
	gdouble  dist;

	//below functions return pointers to actual vertices (whatever
	//the internal surface representation is!
	gts_triangle_vertices (t1, &(v[0]), &(v[1]), &(v[2]));
	gts_triangle_vertices (t2, &(w[0]), &(w[1]), &(w[2]));

	for( i = 0;  i < 3;  i++ )
	{
		equal[i] = 0;
		for( j = 0;  j < 3;  j++ )
		{
		    p1 = &(v[i]->p);     p2 = &(w[j]->p);
		    //CAVEAT - 0.05 works for right now,
		    // should be some epsilon!
		    
		    dist = (gdouble) gts_point_distance2(p1,p2);
		    
		    if( dist < 0.05)
		    {
			equal[i] = 1;

			//FLAG vertices in internal surface representation
			Fi = GTS_OBJECT(v[i])->flags;
			Fj = GTS_OBJECT(w[j])->flags;
			GTS_OBJECT(v[i])->flags |= Fj;
			GTS_OBJECT(w[j])->flags |= Fi;
		    }
		    
		}
	}

	return ( equal[0] && equal[1] && equal[2] );
	//returns 1 if triangles have all vertices matched, 0 otherwise
}


void	flag_initial(GtsObject * o, guint32 *flag)
{
  //guint32 F = *flag;
    o->flags = (guint32)(*flag); //assignment - overrides any previous flag!
}

/**
static void flag_initial_bw(GtsObject * o, guint32 *flag)
{
  //guint32 flag = (guint32)data[0];
    o->flags |= (*flag); //bitwise - does not override previous flag!
}
**/


/* flags triangle with combination (OR, effectively sum) of its vertices' flags */
void	flag_triangle_from_vert(GtsTriangle *t)
{
	GtsVertex *v[3];
	GtsObject *o;
	guint32    v_flag;
	guint	i;

	gts_triangle_vertices (t, &(v[0]), &(v[1]), &(v[2]));
	v_flag = 0;
	for( i = 0;  i < 3;  i++ )
	{
	  o = GTS_OBJECT(v[i]);
	  //printf("%d  ",o->flags);
	  v_flag |= o->flags;
	}
	//printf("    %d \n",v_flag);

	//note - triangle flags don't have to be initialized
	GTS_OBJECT(t)->flags = v_flag;
}


/* watch out, flags edges only if both vertices have the same flag!! */
void	flag_edge_from_vertices(GtsEdge *e)
{
	GtsVertex *v1, *v2;
	guint32    flag1, flag2;
	
	v1 = GTS_SEGMENT(e)->v1;
	v2 = GTS_SEGMENT(e)->v2;
	
	flag1 = GTS_OBJECT(v1)->flags;
	flag2 = GTS_OBJECT(v2)->flags;
	
	if ( flag1 == flag2 ) GTS_OBJECT(e)->flags = flag1;
}


void	build_list_flag(GtsObject *o, Pass_Data *dataptr )
{
  guint32 F1 = dataptr->flag1;
  guint32 F2 = dataptr->flag2;
  
  //printf("\n%d",o->flags);
  if( (o->flags == F1) || (o->flags == F2)  )
  {
	// always use O(1) g_slist_prepend instead of O(n) g_slist_append
	dataptr->list=  g_slist_prepend (dataptr->list, o);
  }
}



GtsSurface *get_surface_from_triangle_flags(GtsSurface *surf,guint32 flag1,guint32 flag2)
{
  GtsSurface *new_surf;
  GSList      *faces = NULL, *f;
  GtsFace    *t;
 

  Pass_Data  data;

  new_surf = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());

			     // build lists of faces

  faces = NULL;  //do not forget to initialize!!

  data.list = faces; //pack data before passing to the func
  data.flag1 = flag1;
  data.flag2 = flag2;
  gts_surface_foreach_face(surf,(GtsFunc) build_list_flag, &data);
  faces = data.list;

  // organize list of faces into a new surface
  f = faces;
  while(f)
  {
	t = GTS_FACE(f->data);
	gts_surface_add_face(new_surf,t);
	f = f->next;
  }

  //{
  //   guint	 n, n1;
  //   n = gts_surface_face_number(surf);
  //   n1 = gts_surface_face_number(new_surf);
  //   printf("\nnumber of faces old %d new %d\n",n,n1); fflush(stdout);
  //} 
  return new_surf;
}


GSList *get_edgelist_from_edge_flags(GtsSurface *surf,guint32 flag1,guint32 flag2)
{  
  GSList      *edges;
 
  Pass_Data  data;
 
  edges = NULL;  //do not forget to initialize!!

  data.list = edges; //pack data before passing to the func
  data.flag1 = flag1;
  data.flag2 = flag2;
  gts_surface_foreach_edge(surf,(GtsFunc) build_list_flag, &data);
  edges = data.list;

  return   edges;
}

GSList *get_vertexlist_from_vertex_flags(GtsSurface *surf,guint32 flag1,guint32 flag2)
{  
  GSList      *vertices;
 
  Pass_Data  data;
 
  vertices = NULL;  //do not forget to initialize!!

  data.list = vertices; //pack data before passing to the func
  data.flag1 = flag1;
  data.flag2 = flag2;
  gts_surface_foreach_vertex(surf,(GtsFunc) build_list_flag, &data);
  vertices = data.list;

  return   vertices;
}


GtsFace *face_from_vertices_safe(GtsVertex *v0, GtsVertex *v1, GtsVertex *v2)
{
   GtsFace *f;
   GtsEdge  *e0,*e1,*e2;
   gdouble  dist0, dist1, dist2;
    
   f = NULL;
   dist0 = gts_point_distance2(GTS_POINT(v0), GTS_POINT(v1));
   if(dist0)
   {
      dist1 = gts_point_distance2(GTS_POINT(v1), GTS_POINT(v2));
      if( dist1 )
      {      
        dist2 = gts_point_distance2(GTS_POINT(v2), GTS_POINT(v0));
      
        if( dist2 )
	{
	  e0 = gts_edge_new(gts_edge_class(),v0,v1);
          e1 = gts_edge_new(gts_edge_class(),v1,v2);
          e2 = gts_edge_new(gts_edge_class(),v2,v0);
   
          f = gts_face_new(gts_face_class(),e0,e1,e2);		
	}
      }
   }      
   
   return f;
}
 

GtsFace *face_from_vertices(GtsVertex *v0, GtsVertex *v1, GtsVertex *v2)
{
   GtsFace *f;
   GtsEdge  *e0,*e1,*e2;
    
   e0 = gts_edge_new(gts_edge_class(),v0,v1);
   e1 = gts_edge_new(gts_edge_class(),v1,v2);
   e2 = gts_edge_new(gts_edge_class(),v2,v0);

   f = gts_face_new(gts_face_class(),e0,e1,e2);
   
   return f;
}
 
   

void  split_surface_from_triangle_flags(
      GtsSurface    *surf,
      guint32       flag1,
      guint32       flag2,
 //     guint32       flag_triple,
      gdouble       iso_fluid,
      GtsSurface    **surf1,
      GtsSurface    **surf2)
{
  GtsSurface *new_surf, *new_surf1, *new_surf2;
  GSList                *faces = NULL, *f;
  GtsFace               *new_face1, *new_face2, *new_face3;
  GtsFace               *ft;
  GtsTriangle           *t;
  Pass_Data             data;
  
  GtsVertex   *v[3], *v_new[3];
  GtsEdge     *e[3];
  GtsObject   *o;
  guint32     v_flag, vertex_flag[3];
  guint	      i, i_next, edge_split[3];
  gdouble     x,y,z;
  gdouble     v1, v2, *pv1, *pv2,rho, fluid_value, *pfluid_value;
  guint32     flag_both = flag1 | flag2;
  guint       add1, add2, add_both;
   
  new_surf1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());
			     
  new_surf2 = gts_surface_new (gts_surface_class (),gts_face_class (),
			     gts_edge_class (),gts_vertex_class ());

  // build lists of faces
  faces = NULL;  //initialize!!

  //pack data before passing to the function
  data.list = faces; 
  data.flag1 = flag1 | flag2;
  data.flag2 = flag2 | flag1;
  gts_surface_foreach_face(surf,(GtsFunc) build_list_flag, &data);
  faces = data.list;

  // organize list of faces into two new surfaces
  f = faces;
  
  /* we assume that our surface (list of faces) is a mixture of triangles whose vertices
     are either marked flag1 or flag2, as it usually happens at boundary of two different
     surfaces */
     
  while(f)
  {
	t = GTS_TRIANGLE(f->data);		
	gts_triangle_vertices_edges (t, NULL, &(v[0]), &(v[1]), &(v[2]),&(e[0]), &(e[1]), &(e[2]) );
	v_flag = 0;
	for( i = 0;  i < 3;  i++ )
	{
	  o = GTS_OBJECT(v[i]);
	  vertex_flag[i] = o->flags;
	  v_flag |= o->flags;
	}
		
        /* along an edge whose vertices are of different types (flags) 
	   we need to
	   add a vertex (whose flag will be both - flag1 | flag2 ) */		
	for( i = 0;  i < 3;  i++ )
	{
	   edge_split[i] = 0;
	   
	   i_next = (i+1)%3;
	   
	   if( vertex_flag[i] != vertex_flag[ i_next ] ) 
	   {	       	        	       
	       pv1 = (gdouble *)(GTS_OBJECT(v[i])->reserved);
	       pv2 = (gdouble *)(GTS_OBJECT(v[i_next])->reserved);
	       	       
	       v1 = *pv1 - iso_fluid;
	       v2 = *pv2 - iso_fluid;
	      
	       if( fabs(v1) < FLT_EPSILON )
	       {
	           GTS_OBJECT(v[i])->flags = flag1 | flag2;
		   vertex_flag[i] = flag1 | flag2;
	       }	       
	       else if ( fabs(v2) < FLT_EPSILON )
	       {
	           GTS_OBJECT(v[i_next])->flags = flag1 | flag2;
		   vertex_flag[i_next] = flag1 | flag2;	       
	       }
	       else
	       {
	          edge_split[i] = 1;
		  rho = v1/(v1-v2);
		  //printf("\nrho %g",rho);
		  x = GTS_POINT(v[i])->x  + rho*(GTS_POINT(v[i_next])->x - GTS_POINT(v[i])->x);
		  y = GTS_POINT(v[i])->y  + rho*(GTS_POINT(v[i_next])->y - GTS_POINT(v[i])->y);
		  z = GTS_POINT(v[i])->z  + rho*(GTS_POINT(v[i_next])->z - GTS_POINT(v[i])->z);

		  v_new[i] = gts_vertex_new (gts_vertex_class(),x,y,z);

		  GTS_OBJECT(v_new[i])->flags = flag1 | flag2;
                 //	       printf("\nflagged %d",GTS_OBJECT(v_new[i])->flags);

		  pfluid_value = (gdouble *)malloc(sizeof(double));
		  *pfluid_value = iso_fluid;
		  GTS_OBJECT(v_new[i])->reserved = pfluid_value;					      

		  //printf("\n%g %g %g split val %g fld_value %g",
		  //      GTS_POINT(v[i])->x,GTS_POINT(v[i])->y,GTS_POINT(v[i])->z,*pv1,fluid_value); 	
	      }
	   }	
	}		
	
	//printf("\nedge_split %d - %d %d %d",edge_split[0]+edge_split[1] + edge_split[2],
	//                       vertex_flag[0], vertex_flag[1], vertex_flag[2]);
	
	/* Now split triangle into 2 or 3 new ones, and add them  to the correct surface */	
	if(edge_split[0] && edge_split[1] )
	{
	    new_face1 = face_from_vertices(v_new[0],v[1],v_new[1]);
	    new_face2 = face_from_vertices(v[0],v_new[0],v_new[1]);
	    new_face3 = face_from_vertices(v_new[1],v[2],v[0]);	    	    
	    
	    if(vertex_flag[1] == flag1 )
	    {
	      gts_surface_add_face(new_surf1,new_face1);
	      
	      gts_surface_add_face(new_surf2,new_face2);
	      gts_surface_add_face(new_surf2,new_face3);
	      
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);
	      
	      gts_surface_add_face(new_surf1,new_face2);
	      gts_surface_add_face(new_surf1,new_face3);
	    }	    	      
	}
	else if(edge_split[1] && edge_split[2] )
	{
	    new_face1 = face_from_vertices(v_new[1],v[2],v_new[2]);
	    new_face2 = face_from_vertices(v[0],v[1],v_new[1]);
	    new_face3 = face_from_vertices(v_new[1],v_new[2],v[0]);

	    if(vertex_flag[2] == flag1  )
	   {
	      gts_surface_add_face(new_surf1,new_face1);
	      
	      gts_surface_add_face(new_surf2,new_face2);
	      gts_surface_add_face(new_surf2,new_face3);
	      
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);
	      
	      gts_surface_add_face(new_surf1,new_face2);
	      gts_surface_add_face(new_surf1,new_face3);
	    }	    	      
	}
	else if(edge_split[2] && edge_split[0] )
	{
	    new_face1 = face_from_vertices(v_new[2],v[0],v_new[0]);
	    new_face2 = face_from_vertices(v_new[0],v[1],v[2]);
	    new_face3 = face_from_vertices(v[2],v_new[2],v_new[0]);
	    
	    if(vertex_flag[0] == flag1 )
	    {
	      gts_surface_add_face(new_surf1,new_face1);
	      
	      gts_surface_add_face(new_surf2,new_face2);
	      gts_surface_add_face(new_surf2,new_face3);
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);
	      
	      gts_surface_add_face(new_surf1,new_face2);
	      gts_surface_add_face(new_surf1,new_face3);
	    }	    	      
	}
	else if(edge_split[0])
	{
	    new_face1 = face_from_vertices(v_new[0],v[1],v[2]);
	    new_face2 = face_from_vertices(v[0],v_new[0],v[2]);    	    
	    
	    if(vertex_flag[1] == flag1 )
	    {
	      gts_surface_add_face(new_surf1,new_face1);	      
	      gts_surface_add_face(new_surf2,new_face2);
	      
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);	      
	      gts_surface_add_face(new_surf1,new_face2);
	    }	    	      	
	}
	else if(edge_split[1])
	{
	    new_face1 = face_from_vertices(v_new[1],v[2],v[0]);
	    new_face2 = face_from_vertices(v[0],v[1],v_new[1]);

	    if(vertex_flag[2] == flag1  )
	   {
	      gts_surface_add_face(new_surf1,new_face1);	      
	      gts_surface_add_face(new_surf2,new_face2);	      
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);	      
	      gts_surface_add_face(new_surf1,new_face2);
	    }	    	      	
	}
	else if(edge_split[2])
	{
	    new_face1 = face_from_vertices(v_new[2],v[0],v[1]);
	    new_face2 = face_from_vertices(v_new[2],v[1],v[2]);
	    
	    if(vertex_flag[0] == flag1 )
	    {
	      gts_surface_add_face(new_surf1,new_face1);	      
	      gts_surface_add_face(new_surf2,new_face2);
	    }  
	    else
	    {
	      gts_surface_add_face(new_surf2,new_face1);	      
	      gts_surface_add_face(new_surf1,new_face2);
	    }	    	      
	}
	else
	{
	    new_face1 = face_from_vertices(v[0],v[1],v[2]);
	    
	    add1 =( (vertex_flag[0] == flag1) ||
	           (vertex_flag[1] == flag1) ||
		   (vertex_flag[2] == flag1) );
		   
	    add2 =( (vertex_flag[0] == flag2) ||
	            (vertex_flag[1] == flag2) ||
		    (vertex_flag[2] == flag2) );		
            		
	    if(add1) 
            {
	      gts_surface_add_face(new_surf1,new_face1);
	    }
	    else if(add2)
	    {  
	      gts_surface_add_face(new_surf2,new_face1);
	    }	    
	}		
		
	
	f = f->next;
  }

 /* {
     guint	 n, n1,n2;
     n = gts_surface_face_number(surf);
     n1 = gts_surface_face_number(new_surf1);
     n2 = gts_surface_face_number(new_surf2);
     printf("\nnumber of faces old %d new1 %d new2 %d\n",n,n1,n2); fflush(stdout);
  }*/ 
  
  *surf1 = new_surf1;
  *surf2 = new_surf2;
}




//Does not free pointer pack!!
void destroy_surface_pack(SurfacePack *pack)
{				       
   /* destroy surfaces we don't need*/
  gts_object_destroy (GTS_OBJECT (pack->surf1));
  gts_object_destroy (GTS_OBJECT (pack->surf2));
  gts_object_destroy (GTS_OBJECT (pack->surfg));
  gts_object_destroy (GTS_OBJECT (pack->surf12));
  gts_object_destroy (GTS_OBJECT (pack->surfg1));
  gts_object_destroy (GTS_OBJECT (pack->surf1g));
  gts_object_destroy (GTS_OBJECT (pack->surfg2));
  gts_object_destroy (GTS_OBJECT (pack->surf_1only));
  gts_object_destroy (GTS_OBJECT (pack->surf_2only));
  gts_object_destroy (GTS_OBJECT (pack->surf_gonly));
  
}

void print_areas(SurfacePack *pack,gchar *base,gboolean web_plot, gint blob_id)
{
  gchar fname[256];
  FILE  *fp;
  
  gdouble  area_12, area_g1, area_g2, area_1only, area_2only, area_gonly;
  gdouble  area_1, vol_1;
  
 /** Compute areas of all surfaces **/  
  area_12 = gts_surface_area (pack->surf12);
  area_g1 = gts_surface_area (pack->surfg1);
  area_g2 = gts_surface_area (pack->surfg2);

  area_1only = gts_surface_area (pack->surf_1only);
  area_2only = gts_surface_area (pack->surf_2only);
  area_gonly = gts_surface_area (pack->surf_gonly);

   
  sprintf(fname,"%s_areas",base);
  fp = fopen(fname,"a");
  if( web_plot) fprintf(fp,"\nBlob %d",blob_id);

  fprintf(fp,"\nAreas S12 %g S1G %g S2G %g S1only %g S2only %g SGonly %g",
                  area_12, area_g1, area_g2,area_1only,area_2only,area_gonly);
  //Area 1 presumed closed (not checked), so output total area & volume
  area_1 =  area_12 + area_g1 + area_1only;
  vol_1 = gts_surface_volume(pack->surf1);
  fprintf(fp,"\nTotal_S1_area %g  volume %g ratio %g",area_1,
	                                                 vol_1,area_1/vol_1);
  fclose(fp);
}




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
    gdouble   a, b, c, rho, one_minus_rho, f0, f1, f2, f_1;
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
          f0 = f[ind];           f_1 = f[ind+shift_z];
	  
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
          f0 = f[ind+shift_z];   f_1 = f[ind];
	  
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
          f0 = f[ind];      f_1 = f[ind+shift_x];
	  
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
          f0 = f[ind+shift_x];   f_1 = f[ind];
	  
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
          f0 = f[ind];    f_1 = f[ind+shift_y];
	  
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
          f0 = f[ind+shift_y];  f_1 = f[ind];
	  
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
    if ( vplus1 < 0. && vplus2 < 0. ) 
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
	   // fluid_value = 0;
	  }  
    }
    else if ( vplus1 < 0.)
    {
	  b = f1 - a;
	  fluid_value = a - b * rho;
	  if( fabs(fluid_value) < TOLERANCE ) fluid_value = 0;
	  printf("\nf %g f0 %g f1 %g", fluid_value,f0,f1);
	  if(( ( fluid_value > f0) && (f0 < f1) ) || (( fluid_value > f0 ) && (f1 > f0) ) ) printf(" test");
          if( (( fluid_value > 0) && (f0 < 0)) || (( fluid_value < 0) && (f0 > 0)) ) 
	  {
	    printf(" test1");
	   // fluid_value = 0;
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




/* 
   vertex_interpolate_value()
   
   This function is designed to work on vertices of an isosurface and interpolate some
   external data on those vertices. Vertices are also classified according to whether the
   interpolated value is below or above iso_interpolate_data.
   interp_value < pdata->iso_interpolate_data --> flag_below
   interp_value >= pdata->iso_interpolate_data --> flag_above
   
   Note: this interpolation is correct only for the vertices located on edges between grid points
         (as is the case for isosurface vertices).
         For points that are strictly inside cubes delineated by grid points, general trilinear interpolation
	 has to be applied.   
*/
void vertex_interpolate_value(GtsVertex *v, Data_Interpolate_Value *pdata)
{
    gint    i, j, k, i1, j1, k1, ind, ind1;
    gint    nx = (pdata->g)->nx, nxy = ((pdata->g)->nx) * ((pdata->g)->ny);    
    
    gdouble   x,y,z;
    gdouble   *d, *f;
    gdouble   d0, d1,f0,f1;
    gdouble   interp_value, *pinterp_value, test;
    gdouble   rho;
    GtsCartesianGrid *g;
    guint32  flag;  
    
    g_return_if_fail (pdata->data_to_interpolate != NULL);
    g_return_if_fail (pdata->data != NULL);
    
    /* shortcuts */
    d = pdata->data;
    f = pdata->data_to_interpolate;
    g = pdata->g;
       
    x =  GTS_POINT(v)->x;
    y =  GTS_POINT(v)->y;
    z =  GTS_POINT(v)->z;

    /* Locate edge that the vertex lies on */
    i =  floor(x - g->x);
    j =  floor(y - g->y);
    k =  floor(z - g->z);
    
    i1 = ceil(x - g->x);
    j1 = ceil(y - g->y);
    k1 = ceil(z - g->z);   
   
    /* find packed indices of edge endpoints in 3d array */
    ind =  k *nxy + j *nx + i; 
    ind1 = k1*nxy + j1*nx + i1;
    
    d0 = d[ind] - pdata->iso;
    f0 = f[ind];
    
    if( (i1 <= (g->nx - 1)) && (j1 <= (g->ny - 1)) && (k1 <= (g->nz - 1)) )
    {
      d1 = d[ind1] - pdata->iso;
      f1 = f[ind1];
      
      if(d1*d0 > 0) printf("\nBummer in vertex_interpolate_value(), vertex might not be on the edge");
      
       /* Interpolation */
      rho = d0/(d0 - d1);
    
      interp_value = f0 + rho*(f1 - f0);

      
      if( fabs( interp_value - pdata->iso_interpolate_data ) < TOLERANCE )
      { 
          interp_value = pdata->iso_interpolate_data;  
      }
    }
    else
    {  
      rho = 0;
      d1 = d0; f1 = f0;
      interp_value = f0;
    }  
       
    pinterp_value = (gdouble *)malloc(sizeof(gdouble));
                      
    test = interp_value - pdata->iso_interpolate_data;
    if ( test < 0)
      flag = pdata->flag_below;
    else //if(test > 0)
      flag = pdata->flag_above;
    //else
      //flag = pdata->flag_above | pdata->flag_below;

    GTS_OBJECT(v)->flags = flag;
    *pinterp_value = interp_value;
    GTS_OBJECT(v)->reserved = pinterp_value;        
}


/* 
   vertex_interpolate_value_color()
   
   This function is designed to work on vertices of an isosurface and interpolate some
   external data on those vertices. Vertices are also classified according to whether the
   interpolated value is below or above iso_interpolate_data.
   interp_value < pdata->iso_interpolate_data --> flag_below
   interp_value >= pdata->iso_interpolate_data --> flag_above
   
   This function also needs to have "data_to_interpolate_color" provided. This data
   is interpolated as well, and will be used to color vertices.
   
   Note: this interpolation is correct only for the vertices located on edges between grid points
         (as is the case for isosurface vertices).
         For points that are strictly inside cubes delineated by grid points, general trilinear interpolation
	 has to be applied.   
*/
void vertex_interpolate_value_color(GtsVertex *v, Data_Interpolate_Value *pdata)
{
    gint    i, j, k, i1, j1, k1, ind, ind1;
    gint    nx = (pdata->g)->nx, nxy = ((pdata->g)->nx) * ((pdata->g)->ny);    
    
    gdouble   x,y,z;
    gdouble   *d, *f, *f_color;
    gdouble   d0, d1,f0,f1, f_color0, f_color1;
    gdouble   interp_value, interp_value_color, *pinterp_value,  test;
    gdouble   rho;
    GtsCartesianGrid *g;
    guint32  flag;  
    
    g_return_if_fail (pdata->data_to_interpolate != NULL);
    g_return_if_fail (pdata->data_to_interpolate_color != NULL);
    g_return_if_fail (pdata->data != NULL);
    
    /* shortcuts */
    d = pdata->data;
    f = pdata->data_to_interpolate;
    f_color = pdata->data_to_interpolate_color;
    g = pdata->g;
       
    x =  GTS_POINT(v)->x;
    y =  GTS_POINT(v)->y;
    z =  GTS_POINT(v)->z;

    /* Locate edge that the vertex lies on */
    i =  floor(x - g->x);
    j =  floor(y - g->y);
    k =  floor(z - g->z);
    
    i1 = ceil(x - g->x);
    j1 = ceil(y - g->y);
    k1 = ceil(z - g->z);   
   
    /* find packed indices of edge endpoints in 3d array */
    ind =  k *nxy + j *nx + i; 
    ind1 = k1*nxy + j1*nx + i1;
    
    d0 = d[ind] - pdata->iso;
    f0 = f[ind];
    f_color0 = f_color[ind];
    
    if( (i1 <= (g->nx - 1)) && (j1 <= (g->ny - 1)) && (k1 <= (g->nz - 1)) )
    {
      d1 = d[ind1] - pdata->iso;
      f1 = f[ind1];
      f_color1 = f_color[ind1];
      
      if(d1*d0 > 0) printf("\nBummer in vertex_interpolate_value(), vertex might not be on the edge");
      
       /* Interpolation */
      rho = d0/(d0 - d1);
    
      interp_value = f0 + rho*(f1 - f0);
      interp_value_color = f_color0 + rho*(f_color1 - f_color0);

      
      if( fabs( interp_value - pdata->iso_interpolate_data ) < TOLERANCE )
      { 
          interp_value = pdata->iso_interpolate_data;  
      }
    }
    else
    {  
      rho = 0;
      d1 = d0; 
      f1 = f0; 
      f_color1 = f_color0;
      interp_value = f0;
      interp_value_color = f_color0;
    }  
       
    if(pdata->use_color_minmax_input == 0)
    { 
      if(interp_value_color >  pdata->color_max) pdata->color_max = interp_value_color;
      if(interp_value_color <  pdata->color_min) pdata->color_min = interp_value_color;  
    }
       
    pinterp_value = (gdouble *)malloc(2*sizeof(gdouble));
                      
    test = interp_value - pdata->iso_interpolate_data;
    if ( test < 0)
      flag = pdata->flag_below;
    else //if(test > 0)
      flag = pdata->flag_above;
    //else
      //flag = pdata->flag_above | pdata->flag_below;

    GTS_OBJECT(v)->flags = flag;
    pinterp_value[0] = interp_value;
    pinterp_value[1] = interp_value_color;
    GTS_OBJECT(v)->reserved = pinterp_value;        
}

static FILE *fp_colors;
int    first = 1;

/* Scales interpolated color value for each vertex to from [cmin,cmax] to [0,1] */
void scale_interpolate_value_color(GtsVertex *v, Data_Interpolate_Value *pdata)
{
  
    gdouble cmin = pdata->color_min;
    gdouble cmax = pdata->color_max;
    gdouble c;
    
    gdouble *pinterp_value = (gdouble *)GTS_OBJECT(v)->reserved;
    c = pinterp_value[1];
    
    if(first)
    {
       fp_colors = fopen("interpolated_values","w");
       first = 0;
    }
    
    fprintf(fp_colors,"%g\n",c);
    
    if (c < cmin) c = 0;
    else if (c > cmax) c = 1;    
    else if ( (fabs(cmax - cmin) < TOLERANCE) ) 
       c = 0.0;
    else
       c = (c - cmin) / (cmax - cmin);
       
    pinterp_value[1] = c;      
}


void interpolate_measure_contact_angle_at_triple_line_edges(GtsSurface *surf_split12, 
                                                GtsSurface *surfg,
						GNode *tree_surfg, 
                                                GSList *triple_line,
						GtsCartesianGrid *g)
{  
  /* We want to 
   a) For each edge on the triple line find the incident triangle from S1 + S2 (i.e. surf_split12 surface).  
   b) Each edge in a) find grain surface triangle that contains/intersects it.
   c) Determine angle that triangles in b) and c) enclose for each edge. */
                  
   
    GSList *l;
    GtsVertex *v;
    GtsEdge *e;
    GtsSegment *s;
    GtsBBox *surfg_triangle_bbox = NULL;
    gdouble distance, normal12[3], normalG1[3];
    gdouble  norm12, normG1, tmp, cos_theta, *ptheta;
    gdouble  theta_avg, num_edges, num_skipped;
    GtsTriangle *t12, *tG1; 
  
    GtsSurface *test12, *testG1, *test12_skip, *testG1_skip;
    
    
    test12 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			       
   
    testG1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());	

    test12_skip = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			       
   
    testG1_skip = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
			    			    
    l = triple_line;
  
    theta_avg = num_edges = num_skipped = 0;
    while(l)
    {
       s = GTS_SEGMENT(l->data);
       e = GTS_EDGE(l->data); /*edge on the triple line */
       
       /* find triangle incident with this edge that is part of S1+S2 surface */
       t12 = GTS_TRIANGLE( gts_edge_is_boundary(e,surf_split12) );
        
       /* get that triangle normal */	            
       gts_triangle_normal(t12,&(normal12[0]),&(normal12[1]),&(normal12[2]));
       
       /* We need to find grain surface triangle that the triple line edge lies in */
       /* We'll use midpoint of the edge for search */
       v = gts_segment_midvertex(s,gts_vertex_class());

       /* if i don't look at the midvertex only, I could in some cases find 2 triangles that are
         cut by the edge */
       distance = gts_bb_tree_point_distance(tree_surfg,GTS_POINT(v),
                                             (GtsBBoxDistFunc)gts_point_triangle_distance2,
					     &surfg_triangle_bbox);
       tG1 = GTS_TRIANGLE(surfg_triangle_bbox->bounded);       
       
       gts_triangle_normal(tG1,&(normalG1[0]),&(normalG1[1]),&(normalG1[2]));

       cos_theta = normal12[0] * normalG1[0] + normal12[1] * normalG1[1] 
        	 + normal12[2] * normalG1[2];
	//normalG1 points outside grain phase, normal12 points outside fluid 1 phase
	// angle(nG1,n12) is (180 - angle_we_need)
       cos_theta*=  -1.0; 
       norm12 = sqrt( normal12[0] * normal12[0] + normal12[1] * normal12[1] 
        	 + normal12[2] * normal12[2]);
       normG1 = sqrt( normalG1[0] * normalG1[0] + normalG1[1] * normalG1[1] 
        	 + normalG1[2] * normalG1[2]);

       tmp = norm12*normG1; 
       if( tmp > 0 )
       {
           cos_theta/= tmp;
           ptheta = (gdouble *)MALLOC(sizeof(gdouble));
           *ptheta = acos(cos_theta) * 180.0 / G_PI;
           //printf("\ntheta %g",*ptheta); 
      
           //use reserved pointer to store pointer to theta
           GTS_OBJECT(e)->reserved = ptheta;

           theta_avg += *ptheta;
           num_edges += 1;
	   
	   gts_surface_add_face(test12,GTS_FACE(t12));
	   gts_surface_add_face(testG1,GTS_FACE(tG1));
       }
       else
       {
           num_skipped +=1;
	   
	   if(norm12 == 0) gts_surface_add_face(test12_skip,GTS_FACE(t12));
	   else if(normG1 == 0) gts_surface_add_face(testG1_skip,GTS_FACE(tG1));
       }

       
       l = l->next;
    }
    
     theta_avg/=num_edges;
     
     printf("\ntheta average: %g num edges %g num_skipped %g",theta_avg,num_edges, num_skipped); 

/*
     gts_surface_foreach_face (testG1, (GtsFunc)set_blue_color, NULL);    
     gts_output_surface_oogl_spec(testG1,g,"interp_testG1.list");
      
     gts_surface_foreach_face (test12, (GtsFunc)set_yellow_color, NULL);
     gts_output_surface_oogl_spec(test12,g,"interp_test12.list");
     
       gts_surface_foreach_face (testG1_skip, (GtsFunc)set_blue_color, NULL);    
     gts_output_surface_oogl_spec(testG1_skip,g,"interp_testG1_skip.list");
      
     gts_surface_foreach_face (test12_skip, (GtsFunc)set_yellow_color, NULL);
     gts_output_surface_oogl_spec(test12_skip,g,"interp_test12_skip.list");
*/     
  }  
  
  

void extrapolate_measure_contact_angle_at_triple_line_edges(GtsSurface *surf_splitG1, 
                                                GtsSurface *surf12,
                                                GSList *triple_line,
						GtsCartesianGrid *g)
{  
  /* We want to 
   a) For each edge on the triple line find the incident triangle from S1 + S2 (i.e. surf_split12 surface).  
   b) Each edge in a) find grain surface triangle that contains/intersects it.
   c) Determine angle that triangles in b) and c) enclose for each edge. */
                  
   
    GSList *l, *surf12_edge_vertices, *surf12_edge_line, *min_ptr, *vert;
    GtsVertex *v;
    GtsEdge *e;
    GtsSegment *s;
    GtsBBox *surfg_triangle_bbox = NULL;
    gdouble distance, normal12[3], normalG1[3];
    gdouble  norm12, normG1, tmp, cos_theta, *ptheta;
    gdouble  theta_avg, num_edges, num_skipped;
    GtsTriangle *t12, *tG1;
    guint32 S1=1, S2 = 2;
  
    GtsSurface *test12, *testG1, *test12_skip, *testG1_skip;
    
    gdouble   dist, min_dist;
    GtsFace *f12;
    
    surf12_edge_line = gts_surface_boundary(surf12);
    gts_output_edgelist_oogl_spec(surf12_edge_line,g,"edge_line.list");
    surf12_edge_vertices = gts_vertices_from_segments(surf12_edge_line);
   
    
    
    test12 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			       
   
    testG1 = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
			    
    test12_skip = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());			       
   
    testG1_skip = gts_surface_new (gts_surface_class (),gts_face_class (),
			    gts_edge_class (),gts_vertex_class ());
			    	
			    			    
    l = triple_line;
  
    theta_avg = num_edges = num_skipped = 0;
    while(l)
    {
       s = GTS_SEGMENT(l->data);
       e = GTS_EDGE(l->data); /*edge on the triple line */
  
  
       vert = surf12_edge_vertices;
       min_dist = 10000;
       
       while(vert)
       {
         dist = gts_point_segment_distance2( GTS_POINT(vert->data), s);
	 if( dist < min_dist)
	 {
	    min_dist = dist;
	    min_ptr = vert;
	 }
         vert = vert->next;
       }
       
       printf("\nmin_dist %g",min_dist); fflush(stdout);
       
       /* triple line, and its segments are inheretnly from grain surface
         mind the orientation of this new triangle */
       f12 = face_from_vertices(s->v1,s->v2,min_ptr->data);
       t12 =  GTS_TRIANGLE(f12);
       
       /* get that triangle normal */	            
       gts_triangle_normal(t12,&(normal12[0]),&(normal12[1]),&(normal12[2]));
            
            /* find triangle incident with this edge that is part of splitt SG+S1 surface */
       tG1 = GTS_TRIANGLE( gts_edge_is_boundary(e,surf_splitG1) );
        
       /* get that triangle normal */	            
       gts_triangle_normal(tG1,&(normalG1[0]),&(normalG1[1]),&(normalG1[2]));
  
       
       cos_theta = normal12[0] * normalG1[0] + normal12[1] * normalG1[1] 
        	 + normal12[2] * normalG1[2];
	//normalG1 points outside grain phase, normal12 points outside fluid 1 phase
	// angle(nG1,n12) is (180 - angle_we_need)
       //cos_theta*=  -1.0; 
       norm12 = sqrt( normal12[0] * normal12[0] + normal12[1] * normal12[1] 
        	 + normal12[2] * normal12[2]);
       normG1 = sqrt( normalG1[0] * normalG1[0] + normalG1[1] * normalG1[1] 
        	 + normalG1[2] * normalG1[2]);

       tmp = norm12*normG1; 
       if( tmp > 0 )
       {
           cos_theta/= tmp;
           ptheta = (gdouble *)MALLOC(sizeof(gdouble));
           *ptheta = acos(cos_theta) * 180.0 / G_PI;
           printf("\ntheta %g",*ptheta); 
      
           //use reserved pointer to store pointer to theta
           GTS_OBJECT(e)->reserved = ptheta;

           theta_avg += *ptheta;
           num_edges += 1;
	   
	   gts_surface_add_face(test12,GTS_FACE(t12));
	   gts_surface_add_face(testG1,GTS_FACE(tG1));
       }
       else
       {
           num_skipped +=1;
	   
	   if(norm12 == 0) gts_surface_add_face(test12_skip,GTS_FACE(t12));
	   else if(normG1 == 0) gts_surface_add_face(testG1_skip,GTS_FACE(tG1));
       }

       
       l = l->next;
    }
    
     theta_avg/=num_edges;
     
     printf("\ntheta average: %g num edges %g num_skipped %g",theta_avg,num_edges, num_skipped); 

     gts_surface_foreach_face (testG1, (GtsFunc)set_blue_color, NULL);    
     gts_output_surface_oogl_spec(testG1,g,"testG1.list");
      
     gts_surface_foreach_face (test12, (GtsFunc)set_yellow_color, NULL);
     gts_output_surface_oogl_spec(test12,g,"test12.list");
     
      gts_surface_foreach_face (testG1_skip, (GtsFunc)set_blue_color, NULL);    
     gts_output_surface_oogl_spec(testG1_skip,g,"testG1_skip.list");
      
     gts_surface_foreach_face (test12_skip, (GtsFunc)set_yellow_color, NULL);
     gts_output_surface_oogl_spec(test12_skip,g,"test12_skip.list");

  }  
  

  
