#ifndef INCLUDE_INTERSECT_SURF_H
#define INCLUDE_INTERSECT_SURF_H


typedef struct _Pass_Data {
    guint32 flag1;
    guint32 flag2;
    GSList *list;
} Pass_Data;

/* Structure for storing 3 surfaces surf1, surf2, and surfg and their
   intersections.
*/
typedef struct _SurfacePack {
  GtsSurface * surf1, *surf2, *surfg;
  GtsSurface *surf12, *surf21;
  GtsSurface *surfg1, *surf1g;
  GtsSurface *surfg2,  *surf2g;
  GtsSurface *surf_1only, *surf_2only, *surf_gonly;
} SurfacePack;

#define SPsz sizeof(SurfacePack)

#define Add_surf(surf,g,fp,txt,set_color,base,fname_add) \
  c = '#'; \
  fprintf(fp,"\n{"); \
  fprintf(fp,"%c%s\n",c,txt); \
  gts_surface_foreach_face (surf, (GtsFunc)set_color,NULL); \
  gts_surface_write_oogl(surf,fp); \
  fprintf(fp,"\n}"); \
\
  sprintf(fname,"%s%s",base,fname_add); \
  gts_output_surface_oogl_spec(surf,&g,fname);


#define Add_surf_only(surf,g,fp,txt,set_color) \
  c = '#'; \
  fprintf(fp,"\n{"); \
  fprintf(fp,"%c%s\n",c,txt); \
  gts_surface_foreach_face (surf, (GtsFunc)set_color,NULL); \
  gts_surface_write_oogl(surf,fp); \
  fprintf(fp,"\n}"); \
  
  
SurfacePack *intersect_surfaces_from_data_base(guchar *, guchar *, guchar *,
                      GtsCartesianGrid *, guint32, guint32, guint32, gboolean, 
                                                       gboolean);

inline gint  match_triangle_vertices(GtsTriangle *,GtsTriangle *);
void intersect_triangles_new(GtsBBox *, GtsBBox *);

void flag_initial(GtsObject *, guint32 *);
void flag_triangle_from_vert(GtsTriangle *);
void flag_edge_from_vertices(GtsEdge *e);


void build_list_flag(GtsObject *, Pass_Data *);
GtsSurface *get_surface_from_triangle_flags(GtsSurface *,guint32,guint32);
GSList *get_edgelist_from_edge_flags(GtsSurface *surf,guint32 flag1,guint32 flag2);
GSList *get_vertexlist_from_vertex_flags(GtsSurface *surf,guint32 flag1,guint32 flag2);

void  split_surface_from_triangle_flags(GtsSurface *surf,guint32 flag1,guint32 flag2,
                                        gdouble iso_fluid,GtsSurface **surf1, GtsSurface **surf2);
					
void print_areas(SurfacePack *,gchar *,gboolean,gint);	     
void destroy_surface_pack(SurfacePack *);	   

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
void vertex_extrapolate_fluid_value(GtsVertex *v, Data_Extrapolate_Fluid_Value *pdata);


/* tolerance 1e-4 works rather well with finney pack test 
  (note that dx=0.04 there; 1e-10 is too little*/
#define TOLERANCE 1e-4


typedef struct _Data_Interpolate_Value{    
    gdouble          *data;                      //(packed) 3d array    
    gdouble          *data_to_interpolate;       //(packed) 3d array
    gdouble          *data_to_interpolate_color; //(packed) 3d array
    gdouble          iso;                    //isolevel
    gdouble          iso_interpolate_data;
    GtsCartesianGrid *g;
    guint32          flag_below;       //flag vertices whose interpolated value is below iso_interpolate_data
    guint32          flag_above;       //flag vertices whose interpolated value is above iso_interpolate_data
                                       //vertices on iso_interpolate_data flagged with (flag_below | flag_above)

    gdouble          color_max;
    gdouble          color_min;
    gdouble          use_color_minmax_input;
} Data_Interpolate_Value;

void vertex_interpolate_value(GtsVertex *v, Data_Interpolate_Value *pdata);
void vertex_interpolate_value_color(GtsVertex *v, Data_Interpolate_Value *pdata);

void scale_interpolate_value_color(GtsVertex *v, Data_Interpolate_Value *pdata);

void interpolate_measure_contact_angle_at_triple_line_edges(GtsSurface *surf_split12, 
                                                GtsSurface *surfg,
						GNode *tree_surfg, 
                                                GSList *triple_line,
						GtsCartesianGrid *g);
	
void extrapolate_measure_contact_angle_at_triple_line_edges(GtsSurface *surf_splitg1, 
                                                GtsSurface *surf12,
                                                GSList *triple_line,
						GtsCartesianGrid *g);
												
#endif
