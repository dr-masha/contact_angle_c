
#ifndef INCLUDE_VIZ_UTIL_H
#define INCLUDE_VIZ_UTIL_H


typedef struct _Color_menu_item {
    gchar     clr;            //color number
    gchar     *name;          //color name
    GtsFunc   func;          //function that sets the color  
} Color_menu_item;

/*** Geomview output business ***/							     
void	print_axes_general(FILE *,gdouble,gdouble,gdouble,gdouble,gdouble,
                                                                     gdouble);
void gts_output_surface_oogl_spec(GtsSurface *,GtsCartesianGrid *,gchar *);

void gts_output_edgelist_oogl_spec(GSList *edges,
				   GtsCartesianGrid *g,gchar *fname);
				   
void gts_surface_write_oogl_color_vertices (GtsSurface * s, FILE * fptr);				   


/**** Object color business ****/
void set_no_color(GtsObject *);

GtsColor interpolated_color(GtsObject *);
void set_interpolated_color(GtsObject * o);

GtsColor red_color(GtsObject *);
void set_red_color(GtsObject *);

GtsColor green_color(GtsObject *);
void set_green_color(GtsObject *);

GtsColor blue_color(GtsObject *);
void set_blue_color(GtsObject *);

GtsColor yellow_color(GtsObject *);
void set_yellow_color(GtsObject *);

GtsColor cyan_color(GtsObject *);
void set_cyan_color(GtsObject *);

GtsColor magenta_color(GtsObject *);
void set_magenta_color(GtsObject *);

#endif
