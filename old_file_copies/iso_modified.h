#ifndef INCLUDE_ISO_MODIFIED_H
#define INCLUDE_ISO_MODIFIED_H

typedef enum { LEFT = 0, RIGHT = 1 } Orientation;

typedef struct {
  GtsVertex * v;
  Orientation orientation;
  
  /* this is modification from the original structure*/
  gdouble     fluid_value;  
} OrientedVertexModified;

typedef struct _GtsIsoSliceModified {
  OrientedVertexModified *** vertices;
  guint nx, ny;
} GtsIsoSliceModified;

GtsIsoSliceModified * gts_iso_slice_modified_new (guint nx, guint ny);


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
				   GtsVertexClass * klassG);

void gts_iso_slice_modified_destroy (GtsIsoSliceModified * slice);
void gts_isosurface_slice_modified (GtsIsoSliceModified * slice1,
			            GtsIsoSliceModified * slice2,
			            GtsSurface * surface);
void gts_isosurface_cartesian_modified (GtsSurface * surfaceG,
			       GtsCartesianGrid g,
			       GtsIsoCartesianFunc f,
			       gpointer data,
			       gpointer data_fluid,
			       gdouble iso,
			       gdouble iso_fluid);				    
				   
#endif
