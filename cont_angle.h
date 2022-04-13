#ifndef INCLUDE_CONTACT_ANGLE_H
#define INCLUDE_CONTACT_ANGLE_H


typedef struct _MasterVertex{
    GtsVertex *v1, *v2, *v3;  //pointer to original vertices
    guint32   flag;
} MasterVertex;

typedef struct _Pass_Data_vca{
    GtsSurface *surf_Aonly;
    GtsSurface *surfAB;
    GtsSurface *surfAC;
    guint32    flagAB; //surfAB vertices flag
    guint32    flagAC; //surfAC vertices flag
    GSList     *thetaA; //list of contact angles
    FILE       *dbg_fp;
} Pass_Data_vca;

#define MVsz  sizeof(MasterVertex)

MasterVertex *create_MasterVertex(void);
gint compare_MasterVertex(MasterVertex *, MasterVertex *);
void build_master_list(GtsVertex *, GSList **);

void calc_and_add_normal(GtsTriangle *,gdouble *);
void vertex_contact_angle1(GtsVertex *, Pass_Data_vca *);
void  find_all_vertex_angles(Pass_Data_vca *);
gdouble get_angle_stats(GSList *, FILE *);
void print_all_angles(GSList *, FILE *, FILE *);
gdouble process_contact_angles(SurfacePack *,GtsCartesianGrid *,guint32,guint32,
                 guint32,gboolean,gboolean,guchar *,gboolean,gint);

#endif
