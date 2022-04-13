#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <float.h>

#include "config.h"
#include "gts.h"
#include "intersect_surf.h"
#include "viz_util.h"
#include "cont_angle.h"
#include "mem.h"


//Find all unique vertices and store the pointers to their origins
//in up to 3 surfaces; if all surfaces considered originate from
//the same surface, vertices are inherited and this is not necessary


MasterVertex *create_MasterVertex(void)
{
    MasterVertex *V;

    V = (MasterVertex *)MALLOC(MVsz);
   
    V->v1 = NULL; V->v2 = NULL; V->v3 = NULL;
    V->flag = 0;
    return V;
}


//TO DO - should I consider vertices at some small epsilon distance equal?
gint compare_MasterVertex(MasterVertex *a, MasterVertex *b)
{
    GtsPoint *p1 =&((a->v1)->p), *p2 =&((b->v1)->p);
    
    if( p1->x < p2->x ) return -1;
    else if (p1->x > p2->x ) return 1;
    else
    {
        if( p1->y < p2->y ) return -1;
        else if (p1->y > p2->y ) return 1;
	else
	{
	    if( p1->z < p2->z ) return -1;
            else if (p1->z > p2->z ) return 1;
	    else 
	    {
	      return 0;
	    }
	}
    }
}


void build_master_list(GtsVertex *v, GSList **v_list)
{
  MasterVertex   *new_v, *elem_v;
  GSList         *elem;

  new_v = create_MasterVertex();
  new_v->v1 = v;
  new_v->flag |= GTS_OBJECT(v)->flags;
  
  elem = g_slist_find_custom(*v_list,new_v,(GCompareFunc)compare_MasterVertex);
  
  if( elem == NULL ) 
  {
      *v_list = g_slist_insert_sorted(*v_list,new_v,
                                           (GCompareFunc)compare_MasterVertex);		    
  }
  else
  {  
       FREE(new_v);
       elem_v = elem->data;
       
       //existing list element should store pointer to vertex v
       if( elem_v->v2 == NULL ) elem_v->v2 = v;
       else if ( elem_v->v3 == NULL ) elem_v->v3 = v;
       else printf("\nbuild_master_list() - need more storage space!!");

       elem_v->flag |= GTS_OBJECT(v)->flags; //update flags
       
       //printf("\nsame vertices v %x v1 %x v2 %x",v,elem_v->v1, elem_v->v2);
  }
}


void calc_and_add_normal(GtsTriangle *t,gdouble *normal)
{
    gdouble x,y,z;
 
    gts_triangle_normal(t,&x,&y,&z);
    normal[0]+= x;
    normal[1]+= y;
    normal[2]+= z;
}


/** GtsVertex *v assumed to be from surf_Aonly
*    A and B are thought of as fluids in this context, C as grain phase.
*    Note that surf_Aonly, surfAB and surfAC are assumed to be derived from the
*    same surface surfA so that vertices are inherited from surfA.
*    Because of inheritance, triangle orientations are consistent with surfA as well
*    (most probably all normals point outside phase A region).
**/

void vertex_contact_angle1(GtsVertex *v, Pass_Data_vca *pdata)
{
    GtsSurface *surf_Aonly = pdata->surf_Aonly;
    GtsSurface *surfAB     = pdata->surfAB; 
    GtsSurface *surfAC     = pdata->surfAC;
    guint32    flagAB      = pdata->flagAB;
    guint32    flagAC      = pdata->flagAC;
    FILE       *dbg_fp     = pdata->dbg_fp;
    
    GSList   *tAB = NULL, *tAC = NULL, *W = NULL, *pW;
    guint    nb, nc;
    gdouble  normalAB[3], normalAC[3], normAB, normAC, tmp;
    GtsVertex   *w;
    gdouble     cos_theta, *ptheta;
    gint   test;
    
    //process vertices (on surf_Aonly) that belong to surfAB as well
    if( GTS_OBJECT(v)->flags == flagAB )
    {
      tAB = gts_vertex_faces(v,surfAB,tAB); //surfAB triangles incident to v
      nb = g_slist_length(tAB);
      //printf("\nnb %d",nb); fflush(stdout);

      normalAB[0] = normalAB[1] = normalAB[2] = 0; 
      if( nb != 0 )
      {
          //normalAB will store sum of all triangle normals   
	  g_slist_foreach(tAB,(GFunc)calc_and_add_normal,normalAB);
          //find average of normals
          normalAB[0]/= (gdouble)nb;
          normalAB[1]/= (gdouble)nb;
          normalAB[2]/= (gdouble)nb;
	  //cntAB++;
      }
      
      //printf("\nnormal %g %g %g,",normalAB[0],normalAB[1], normalAB[2]);
      fflush(stdout);
     
      W = gts_vertex_neighbors(v,W,surf_Aonly);
      
      //nc = g_slist_length(W);
      //printf("\nnW %d flags ",nc); fflush(stdout);

      pW = W;
      while( pW != NULL )
      {
           w = pW->data;
	   
	   //we're interested only in vertices on W that are on surfAC
	   //printf("%d ",GTS_OBJECT(w)->flags);
           if( GTS_OBJECT(w)->flags == flagAC ) 
	   {
              tAC = gts_vertex_faces(w,surfAC,tAC);
           }
	   pW = pW->next;
      }
      //printf("\n");

      nc = g_slist_length(tAC);
      //printf("\nnc %d",nc); fflush(stdout);

      normalAC[0] = normalAC[1] = normalAC[2] = 0;
      if( nc != 0 )
      {
          //normalAC will store sum of all triangle normals
          g_slist_foreach(tAC,(GFunc)calc_and_add_normal,normalAC);
          //find average
          normalAC[0]/= (gdouble)nc;
          normalAC[1]/= (gdouble)nc;
          normalAC[2]/= (gdouble)nc;
	  //cntAC++;
      }
      //printf("\nnormal %g %g %g,",normalAC[0],normalAC[1], normalAC[2]);
      //fflush(stdout);

      if( dbg_fp ) //print normals to Geomview file
      {
          fprintf(dbg_fp,"\n{COFF");
	  fprintf(dbg_fp,"\n3 2 0");
	  fprintf(dbg_fp,"\n%g %g %g 0 0 0 1", (v->p).x,(v->p).y,(v->p).z);
	  fprintf(dbg_fp,"\n%g %g %g 1 0 0 1", (v->p).x + normalAB[0],
	                                     (v->p).y + normalAB[1],
					     (v->p).z + normalAB[2]);
	  fprintf(dbg_fp,"\n%g %g %g 0 0 1 1", (v->p).x - normalAC[0],
	                                     (v->p).y - normalAC[1],
					     (v->p).z - normalAC[2]);
	  fprintf(dbg_fp,"\n2 0 1");
	  fprintf(dbg_fp,"\n2 0 2");
	  fprintf(dbg_fp,"\n}");
      }
      
      //calculate angle between normals
      cos_theta = normalAB[0] * normalAC[0] + normalAB[1] * normalAC[1] 
        	 + normalAB[2] * normalAC[2];
	//becuase of inheritance -- see above assumptions
	//need to reverse normalAC to point towards
	//phase A (and effectively be normal to grain (C) surface
      cos_theta*=  -1.0; 
      normAB = sqrt( normalAB[0] * normalAB[0] + normalAB[1] * normalAB[1] 
        	 + normalAB[2] * normalAB[2]);
      normAC = sqrt( normalAC[0] * normalAC[0] + normalAC[1] * normalAC[1] 
        	 + normalAC[2] * normalAC[2]);

      tmp = normAB*normAC;      
       
      if( tmp > 0 ) 
      {
        cos_theta/= tmp;
	
	/* 4/20/2007 put in cos correction */
	if( cos_theta <= -1) cos_theta = -1.0;
	else if( cos_theta >= 1) cos_theta = 1.0;
	
	ptheta = (gdouble *)MALLOC(sizeof(gdouble));
        *ptheta = acos(cos_theta) * 180.0 / G_PI;
        printf(" theta %g",*ptheta); 
      
        //use reserved pointer to store pointer to theta
        GTS_OBJECT(v)->reserved = ptheta;
      }
      else 
      { 
        /* New fix, 4/20/2007 */
         printf("\nBummer - vertex_contact_angle1 nb %d nc %d",nb,nc);
	 GTS_OBJECT(v)->reserved = NULL; 
      }

      
      
#ifdef STORE_IN_LIST     
      //store in the list as well
      pdata->thetaA =  g_slist_prepend(pdata->thetaA,ptheta);
#endif      
    }
}



void  find_all_vertex_angles(Pass_Data_vca *vca)
{
#ifdef STORE_IN_LIST
  gint     nt;
  GSList   *l;
  gdouble  t, mint, maxt, sumt, *tt;
#endif 
  
  gts_surface_foreach_vertex(vca->surf_Aonly,(GtsFunc)vertex_contact_angle1,
                                                               vca);							 	
#ifdef STORE_IN_LIST  
  l = vca->thetaA;
  mint = 10000.0;
  maxt = -10000.0;
  sumt = 0.0;
  while(l)
  {
     tt = l->data;
     t = *tt;
     //printf("\n%g",t);
     if( t < mint ) mint = t;
     if( t > maxt ) maxt = t;
     sumt+= t;
     l = l->next;
  }
  nt = g_slist_length(vca->thetaA);
  printf("\ntheta (%d values) min %g max %g ave %g\n",nt,mint,maxt,sumt/nt);
  
#endif
}



gdouble get_angle_stats(GSList *master_vert_list, FILE *fp)
{ 
    GSList *l;
    MasterVertex *vM;
    GtsVertex *v1, *v2;
    gdouble t1, t2, *pt1, *pt2, scale_t1, scale_t2;
    GtsRange  theta1, theta2, theta12, scaled1, scaled2;
    FILE *fp1 = stdout;
  
    gts_range_init (&theta1);    gts_range_init (&scaled1); 
    gts_range_init (&theta2);    gts_range_init (&scaled2); 
    gts_range_init (&theta12);
    
    l = master_vert_list;
    while (l)
    { 
	vM = l->data;
	v1 = vM->v1;
	v2 = vM->v2;

	if(( v1 != NULL ) && (v2 != NULL) && 
	    (GTS_OBJECT(v1)->reserved != NULL) &&
	    (GTS_OBJECT(v2)->reserved != NULL ) )
	{
	      pt1 = (gdouble *)(GTS_OBJECT(v1)->reserved);
	      pt2 = (gdouble *)(GTS_OBJECT(v2)->reserved);
	      
	      if( (pt1 != NULL) && (pt2 != NULL) )
	      {
		t1 = *pt1; t2 = *pt2;
		
		scale_t1 = (t1 / (t1+t2)) * 180.0;
		scale_t2 = (t2 / (t1+t2)) * 180.0;
		printf("\n%g %g  %g  scaled %g %g",t1,t2,t1+t2,scale_t1,
		                                                 scale_t2);
		gts_range_add_value (&theta1,t1);
		gts_range_add_value (&scaled1,scale_t1);
		gts_range_add_value (&theta2,t2);
		gts_range_add_value (&scaled2,scale_t2);
		gts_range_add_value (&theta12,t1+t2);
	      } 
	}
	l = l->next;
    }
    
    gts_range_update(&theta1);    gts_range_update(&scaled1);
    gts_range_update(&theta2);    gts_range_update(&scaled2);
    gts_range_update(&theta12);
    
    if( fp != NULL) fp1 = fp;
    
    fprintf(fp1,"\nAngle computed at %d vertices.",theta1.n);
    fprintf(fp1,"\ntheta1 ");
    gts_range_print(&theta1, fp1);
    fprintf(fp1,"\ntheta2 ");
    gts_range_print(&theta2, fp1);
    fprintf(fp1,"\ntheta1 + theta2 ");
    gts_range_print(&theta12,fp1);
    fprintf(fp1,"\nscaled1 ");
    gts_range_print(&scaled1, fp1);
    fprintf(fp1,"\nscaled2 ");
    gts_range_print(&scaled2, fp1);
    fprintf(fp1,"\n");
    
    return theta1.mean; //returns mean of the angle theta1
    
}


void print_all_angles(GSList *master_vert_list, FILE *fp1, FILE *fp2)
{ 
    GSList *l;
    MasterVertex *vM;
    GtsVertex *v1, *v2;
    gdouble t1, t2, *pt1, *pt2;
  
    l = master_vert_list;
    while (l)
    { 
	vM = l->data;
	v1 = vM->v1;
	v2 = vM->v2;

	if(( v1 != NULL ) && (v2 != NULL))
	{
	      pt1 = (gdouble *)(GTS_OBJECT(v1)->reserved);
	      pt2 = (gdouble *)(GTS_OBJECT(v2)->reserved);
	      
	      if( (pt1 != NULL) && (pt2 != NULL) )
	      {
		t1 = *pt1; t2 = *pt2;
		fprintf(fp1,"%g\n",t1);
		fprintf(fp2,"%g\n",t2);
	      } 
	}
	l = l->next;
    }    
}


gdouble process_contact_angles(
     SurfacePack      *pack,      //structure containing all surfaces and their intersections
     GtsCartesianGrid *g,
     guint32      S1,             //surf1 bitwise mark
     guint32      S2,             //surf2 bitwise mark
     guint32      SG,             //surfg bitwise mark
     gboolean     dbg_plot,       //dbg_plot==TRUE will enable plot of pointwise normals
     gboolean     print_pointwise, //print_pointwise==TRUE will print contact angles from all individual points
                                   //otherwise only summary stats are recorded (mean,min,max etc.)
     guchar      *base,            //filenames base, used when (dbg_plot || print_pointwise)
     gboolean     web_plot,
     gint         blob_id)         //used when web_plot==TRUE
{
  Pass_Data_vca  vca1, vca2;
  GSList   *master_vert_list;
  
  FILE *fp, *fp1, *fp2;
  gchar fname[256];
  gdouble mean;
  
/***
  * pointer to angle values for applicable vertices will be stored in 'reserved'
  * so we reset the values beforehand
  * DO NOT USE reserved field for anything else!!!!
  ***/
  gts_surface_foreach_vertex(pack->surf_1only,(GtsFunc)gts_object_reset_reserved,NULL);
  gts_surface_foreach_vertex(pack->surf_2only,(GtsFunc)gts_object_reset_reserved,NULL);

  
  //initialize Pass_data_vca structure
  vca1.surf_Aonly = pack->surf_1only;
  vca1.surfAB = pack->surf12;
  vca1.surfAC = pack->surf1g;
  vca1.flagAB = S1 + S2;
  vca1.flagAC = S1 + SG;
  vca1.thetaA = NULL; // if used to store a list, has to be initialized to NULL
  
  if( dbg_plot )
  {
     sprintf(fname,"%s_dbg_normals.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g->x,g->y,g->z,g->x+g->nx,g->y+g->ny,g->z+g->nz);
     vca1.dbg_fp = fp;
  }
  else
  {
    vca1.dbg_fp = NULL;
  }
    
  find_all_vertex_angles(&vca1);
  if( dbg_plot ) fclose(fp);
  
 
  //initialize Pass_data_vca structure
  vca2.surf_Aonly = pack->surf_2only;
  vca2.surfAB = pack->surf21;
  vca2.surfAC = pack->surf2g;
  vca2.flagAB = S2 + S1;
  vca2.flagAC = S2 + SG;
  vca2.thetaA = NULL; 
  vca2.dbg_fp = NULL;
  
  find_all_vertex_angles(&vca2);
  
  //printf("\ntheta2 pass cntAB %d cntAC %d",cntAB,cntAC);
   
  master_vert_list = NULL;
  gts_surface_foreach_vertex(pack->surf_1only,(GtsFunc)build_master_list,
                                                       &master_vert_list); 
  gts_surface_foreach_vertex(pack->surf_2only,(GtsFunc)build_master_list,
                                                       &master_vert_list);
   //gather statistics on angles
  sprintf(fname,"%s_cont_angle",base);
  fp  = fopen(fname,"a");
  if( web_plot) fprintf(fp,"\nBlob %d",blob_id);
  mean = get_angle_stats(master_vert_list,fp);
  fclose(fp);
  
  if( print_pointwise )
  {//output all theta1 and theta2 pointwise values to separate files
      
      sprintf(fname,"%s_%d_theta1",base,blob_id);
      fp1 = fopen(fname,"w");
      sprintf(fname,"%s_%d_theta2",base,blob_id);
      fp2 = fopen(fname,"w");
      print_all_angles(master_vert_list,fp1,fp2);
      fclose(fp1);
      fclose(fp2);
  }
  
  return mean;
}
