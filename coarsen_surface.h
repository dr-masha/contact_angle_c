#ifndef INCLUDE_COARSEN_SURFACE_H
#define INCLUDE_COARSEN_SURFACE_H

typedef enum { NUMBER, COST } StopOptions;
typedef enum { COST_LENGTH, COST_OPTIMIZED, COST_ANGLE } CostOptions;
typedef enum { MIDVERTEX, OPTIMIZED } MidvertexOptions;

void gts_surface_coarsen_top(GtsSurface * s,gboolean progressive,CostOptions cost,
     MidvertexOptions mid,gboolean log_cost,StopOptions stop );

#endif
