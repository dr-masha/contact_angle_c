#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>

#include "config.h"
#include "gts.h"
#include "viz_util.h"
#include "intersect_surf.h"
#include "mem.h"

/*** Geomview output business ***/

void	print_axes_general(
	FILE	*fp,
	gdouble	xs,
	gdouble	ys,
	gdouble	zs,
	gdouble	xe,
	gdouble	ye,
	gdouble	ze)
{
	fprintf(fp,"{ VECT\n");
	fprintf(fp,"3 6 1\n");
	fprintf(fp,"2 2 2\n");
	fprintf(fp,"1 0 0\n");

	fprintf(fp,"%g %g %g   %g %g %g\n",xs, ys, zs, xe, ys, zs);
	fprintf(fp,"%g %g %g   %g %g %g\n",xs, ys, zs, xs, ye, zs);
	fprintf(fp,"%g %g %g   %g %g %g\n",xs, ys, zs, xs, ys, ze);
	fprintf(fp,"\n0  0  0.75  1\n");		// blue
	fprintf(fp,"}\n");
}

void gts_output_surface_oogl_spec(GtsSurface *s,
				  GtsCartesianGrid *g,gchar *fname)
{
  FILE    * fp;

  fp = fopen(fname,"w");
  fprintf(fp,"LIST\n");
  print_axes_general(fp,g->x,g->y,g->z,g->x+g->nx,g->y+g->ny,g->z+g->nz);
  fprintf(fp,"\n{");
  gts_surface_write_oogl(s,fp);
  fprintf(fp,"\n}");
  fclose (fp);
}



void gts_output_edgelist_oogl_spec(GSList *edges,
				   GtsCartesianGrid *g,gchar *fname)
{
  FILE    * fp;
  GSList   *l;
  GtsSegment *s;
  gint  n;
  
  fp = fopen(fname,"w");
  fprintf(fp,"LIST\n");
  print_axes_general(fp,g->x,g->y,g->z,g->x+g->nx,g->y+g->ny,g->z+g->nz);
  fprintf(fp,"\n{");
  
  l = edges;  
  
  while(l)
  {
     s = GTS_SEGMENT(l->data);
     
     if (GTS_OBJECT (s)->klass->color) 
     {
       GtsColor c = (* GTS_OBJECT (s)->klass->color) (GTS_OBJECT (s));
       fprintf (fp, "VECT 1 2 1 2 1 %g %g %g %g %g %g %g %g %g 1.\n",
             GTS_POINT (s->v1)->x, GTS_POINT (s->v1)->y, GTS_POINT (s->v1)->z,
             GTS_POINT (s->v2)->x, GTS_POINT (s->v2)->y, GTS_POINT (s->v2)->z,
             c.r, c.g, c.b);
     }
    else
    {
      fprintf (fp, "VECT 1 2 0 2 0 %g %g %g %g %g %g\n",
             GTS_POINT (s->v1)->x, GTS_POINT (s->v1)->y, GTS_POINT (s->v1)->z,
             GTS_POINT (s->v2)->x, GTS_POINT (s->v2)->y, GTS_POINT (s->v2)->z);
    }
    
    l  = l->next;
  }
  
  fprintf(fp,"\n}");
  fclose (fp);
}

/**** GtsObject color business ****/

GtsColor    red = {0.8,  0,0}, green = {0,0.8,  0},    blue = {  0,0,0.8};
GtsColor yellow = {0.8,0.8,0},  cyan = {0,0.8,0.8}, magenta = {0.8,0,0.8};

static GtsColor interpolated = {0, 0.4, 0};

Color_menu_item color_menu[] = 
{
    {'0',"none",   (GtsFunc)set_no_color},
    {'1',"red",    (GtsFunc)set_red_color},
    {'2',"green",  (GtsFunc)set_green_color},
    {'3',"yellow", (GtsFunc)set_yellow_color},
    {'4',"blue",   (GtsFunc)set_blue_color},
    {'5',"magenta",(GtsFunc)set_magenta_color},
    {'6',"cyan",   (GtsFunc)set_cyan_color}
};

void set_no_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = NULL;
}

void set_interpolated_color(GtsObject * o)
{
  /* this case assumes that there is a color value in [0,1] that reserved points to */
  /* see vertex_interpolate_value_color() */
  
  gdouble *pinterp_value = (gdouble *)GTS_OBJECT(o)->reserved;
  gdouble red = pinterp_value[1];
  
  interpolated.r = red;
  
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = interpolated_color;
}

/** doesn't work for some reason 

void set_color(GtsObject * o, GtsColor *clr)
{
   // overload color definition
   GTS_OBJECT (o)->klass->color = *clr;
  // sizeof(GTS_OBJECT(surf1)->klass->color));
}
**/

GtsColor interpolated_color(GtsObject * object)
{
  return interpolated;
}


GtsColor red_color(GtsObject * object)
{
  return red;
}


void set_red_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = red_color;
}


GtsColor green_color(GtsObject * object)
{
  return green;
}

void set_green_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = green_color;
}

GtsColor blue_color(GtsObject * object)
{
  return blue;
}

void set_blue_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = blue_color;
}

GtsColor yellow_color(GtsObject * object)
{
  return yellow;
}

void set_yellow_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = yellow_color;
}

GtsColor cyan_color(GtsObject * object)
{
  return cyan;
}

void set_cyan_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = cyan_color;
}


GtsColor magenta_color(GtsObject * object)
{
  return magenta;
}

void set_magenta_color(GtsObject * o)
{
   /* overload color definition */
  GTS_OBJECT (o)->klass->color = magenta_color;
}

#define Add_surf1(surf,g,fp,txt,clr,base) \
  c = '#'; \
  fprintf(fp,"\n{"); \
  fprintf(fp,"%c%s\n",c,txt); \
  gts_surface_foreach_face (surf, (GtsFunc)color_menu[clr].func,NULL); \
  gts_surface_write_oogl(surf,fp); \
  fprintf(fp,"\n}"); \
\
  sprintf(fname,"%s_%s.list",base,color_menu[clr].name);\
  gts_output_surface_oogl_spec(surf,&g,fname);

//CLR1, CLR2, CLRG should be 1(red), 2(green) or 4(blue)
void plot_surface_pack(
     SurfacePack      *pack,
     GtsCartesianGrid  *g,
     guint32          CLR1,
     guint32          CLR2,
     guint32          CLRG,
     gchar            *base)
{
     gchar    fname[256], c;
     FILE    *fp;
     
     guint32   clr;
     
     /* Geomview files of (colored) surfaces are output by default.
      * Each file has local axes added so they can be superimposed at will.
      */ 

     //individual surfaces
     gts_surface_foreach_face(pack->surf1,(GtsFunc)(color_menu[CLR1].func),NULL);

     sprintf(fname,"%s_fld1.list",base);
     gts_output_surface_oogl_spec(pack->surf1,g,fname);

     gts_surface_foreach_face(pack->surf2,(GtsFunc)(color_menu[CLR2].func),NULL);
     sprintf(fname,"%s_fld2.list",base);
     gts_output_surface_oogl_spec(pack->surf2,g,fname);

     gts_surface_foreach_face (pack->surfg,(GtsFunc)(color_menu[CLRG].func), NULL);
     sprintf(fname,"%s_grain.list",base);
     gts_output_surface_oogl_spec(pack->surfg,g,fname);

     //intersections file
     sprintf(fname,"%s_intersect.list",base);
     fp = fopen(fname,"w");
     fprintf(fp,"LIST\n");
     print_axes_general(fp,g->x,g->y,g->z,g->x+g->nx,g->y+g->ny,g->z+g->nz);

     //gts_surface_write(surf12,stdout);
     
     clr = CLR1 + CLR2;
     Add_surf1(pack->surf12,(*g),fp,"S1 & S2",clr,base)
     clr = CLR1 + CLRG;
     Add_surf1(pack->surf1g,(*g),fp,"S1 & SG",clr,base)
     clr = CLR2 + CLRG;
     Add_surf1(pack->surfg2,(*g),fp,"SG & S2",clr,base)

     Add_surf1(pack->surf_1only,(*g),fp,"S1 \\ (S2 U SG)",CLR1,base)
     Add_surf1(pack->surf_2only,(*g),fp,"S2 \\ (S1 U SG)",CLR2,base)
     Add_surf1(pack->surf_gonly,(*g),fp,"SG \\ (S1 U S2)",CLRG,base)
     fclose (fp);
}



/* Variation in surface oogl output in order to get colored vertices */

static void write_vertex_oogl_interpolate_color (GtsPoint * p, gpointer * data)
{
  FILE * fp = data[0];

  fprintf (fp, "%g %g %g", p->x, p->y, p->z);
  set_interpolated_color( GTS_OBJECT (p) );
  
 // if (GTS_OBJECT (p)->klass->color) {
    GtsColor c = (* GTS_OBJECT (p)->klass->color) (GTS_OBJECT (p));
    fprintf (fp, " %g %g %g 1.0\n", c.r, c.g, c.b);
 //}
 // else
    fputc ('\n', fp);
  GTS_OBJECT (p)->reserved = GUINT_TO_POINTER ((*((guint *) data[1]))++);
}

static void write_face_oogl (GtsTriangle * t, FILE * fp)
{
  GtsVertex * v1, * v2, * v3;
  gts_triangle_vertices (t, &v1, &v2, &v3);
  fprintf (fp, "3 %u %u %u",
	   GPOINTER_TO_UINT (GTS_OBJECT (v1)->reserved),
	   GPOINTER_TO_UINT (GTS_OBJECT (v2)->reserved),
	   GPOINTER_TO_UINT (GTS_OBJECT (v3)->reserved));
  if (GTS_OBJECT (t)->klass->color) {
    GtsColor c = (* GTS_OBJECT (t)->klass->color) (GTS_OBJECT (t));
    fprintf (fp, " %g %g %g\n", c.r, c.g, c.b);
  }
  else
    fputc ('\n', fp);
}

/**
 * gts_surface_write_oogl:
 * @s: a #GtsSurface.
 * @fptr: a file pointer.
 * 
 * Writes in the file @fptr an OOGL (Geomview) representation of @s.
 */
void gts_surface_write_oogl_color_vertices (GtsSurface * s, FILE * fptr)
{
  guint n = 0;
  gpointer data[2];
  GtsSurfaceStats stats;

  g_return_if_fail (s != NULL);
  g_return_if_fail (fptr != NULL);

  data[0] = fptr;
  data[1] = &n;

  gts_surface_stats (s, &stats);
  fputs ("COFF ", fptr);
  fprintf (fptr, "%u %u %u\n", 
	   stats.edges_per_vertex.n, 
	   stats.n_faces,
	   stats.faces_per_edge.n);
  gts_surface_foreach_vertex (s, (GtsFunc) write_vertex_oogl_interpolate_color, data);
  gts_surface_foreach_face (s, (GtsFunc) write_face_oogl, fptr);
  gts_surface_foreach_vertex (s, (GtsFunc) gts_object_reset_reserved, NULL);
}
