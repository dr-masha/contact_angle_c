#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include <glib.h>
#include "config.h"
#include "gts.h"


#include "coarsen_surface.h"

/* function for coarsening a surface
 * adapted from coarsen.c example in GTS 
 * GTS - Library for the manipulation of triangulated surfaces
 * Copyright (C) 1999 Stéphane Popinet
 */




static gboolean stop_number_verbose (gdouble cost, guint number, guint * min)
{
  static guint nmax = 0, nold = 0;
  static GTimer * timer = NULL, * total_timer = NULL;

  g_return_val_if_fail (min != NULL, TRUE);

  if (timer == NULL) {
    nmax = nold = number;
    timer = g_timer_new ();
    total_timer = g_timer_new ();
    g_timer_start (total_timer);
  }

  if (number != nold && number % 121 == 0 &&
      number < nmax && nmax > *min) {
    gdouble total_elapsed = g_timer_elapsed (total_timer, NULL);
    gdouble remaining;
    gdouble hours, mins, secs;
    gdouble hours1, mins1, secs1;

    g_timer_stop (timer);

    hours = floor (total_elapsed/3600.);
    mins = floor ((total_elapsed - 3600.*hours)/60.);
    secs = floor (total_elapsed - 3600.*hours - 60.*mins);

    remaining = total_elapsed*((nmax - *min)/(gdouble) (nmax - number) - 1.);
    hours1 = floor (remaining/3600.);
    mins1 = floor ((remaining - 3600.*hours1)/60.);
    secs1 = floor (remaining - 3600.*hours1 - 60.*mins1);

    fprintf (stderr, 
	     "\rEdges: %10u %3.0f%% %6.0f edges/s "
	     "Elapsed: %02.0f:%02.0f:%02.0f "
	     "Remaining: %02.0f:%02.0f:%02.0f ",
	     number, 
	     100.*(nmax - number)/(nmax - *min),
	     (nold - number)/g_timer_elapsed (timer, NULL),
	     hours, mins, secs,
	     hours1, mins1, secs1);
    fflush (stderr);

    nold = number;
    g_timer_start (timer);
  }
  if (number < *min) {
    g_timer_destroy (timer);
    g_timer_destroy (total_timer);
    return TRUE;
  }
  return FALSE;
}

static gboolean stop_cost_verbose (gdouble cost, guint number, gdouble * max)
{
  g_return_val_if_fail (max != NULL, TRUE);

  if (number % 121 == 0) {
    fprintf (stderr, "\rEdges: %10u Cost: %10g ", number, cost);
    fflush (stderr);
  }
  if (cost > *max)
    return TRUE;
  return FALSE;
}

static gboolean stop_log_cost (gdouble cost, guint number)
{
  fprintf (stderr, "%d %g\n", number, cost);
  return FALSE;
}

static gdouble cost_angle (GtsEdge * e)
{
  if (e->triangles && e->triangles->next)
    return fabs (gts_triangles_angle (e->triangles->data,
				      e->triangles->next->data));
  return 100000.0;
}

/* coarsen - produce a coarsened version of the input */
void gts_surface_coarsen_top(
     GtsSurface * s,
     gboolean progressive,
     CostOptions cost,
     MidvertexOptions mid,
     gboolean log_cost,
     StopOptions stop )
{
  GtsPSurface * ps = NULL;
  gboolean verbose = FALSE;
  guint number = 3*gts_surface_edge_number(s)/4; /* reduce number of edges to 3/4
                                  of the original by default */
  gdouble cmax = 0.05; 
  GtsKeyFunc cost_func = NULL;
  GtsCoarsenFunc coarsen_func = NULL;
  GtsStopFunc stop_func = NULL;
  gpointer stop_data = NULL;
  
  gdouble fold = M_PI/180.;
  GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };
  gpointer coarsen_data = NULL, cost_data = NULL;

  /* if verbose on print stats */
  if (verbose) {
    gts_surface_print_stats (s, stderr);
    fprintf (stderr, "# volume: %g area: %g\n", 
	     gts_surface_volume (s), gts_surface_area (s));
  }

  /* select the right coarsening process */
  switch (cost) {
  case COST_OPTIMIZED: 
    cost_func = (GtsKeyFunc) gts_volume_optimized_cost; 
    cost_data = &params;
    break;
  case COST_LENGTH:
    cost_func = NULL; break;
  case COST_ANGLE:
    cost_func = (GtsKeyFunc) cost_angle; break;
  default:
    g_assert_not_reached ();
  }
  switch (mid) {
  case MIDVERTEX:
    coarsen_func = NULL; 
    break;
  case OPTIMIZED:
    coarsen_func = (GtsCoarsenFunc) gts_volume_optimized_vertex; 
    coarsen_data = &params;
    break;
  default:
    g_assert_not_reached ();
  }
  if (log_cost)
    stop_func = (GtsStopFunc) stop_log_cost;
  else {
    switch (stop) {
    case NUMBER:
      if (verbose)
	stop_func = (GtsStopFunc) stop_number_verbose;
      else
	stop_func = (GtsStopFunc) gts_coarsen_stop_number; 
      stop_data = &number;
      break;
    case COST:
      if (verbose)
	stop_func = (GtsStopFunc) stop_cost_verbose;
      else
	stop_func = (GtsStopFunc) gts_coarsen_stop_cost; 
      stop_data = &cmax;
      break;
    default:
      g_assert_not_reached ();
    }
  }
  if (progressive)
    ps = gts_psurface_new (gts_psurface_class (),
			   s, gts_split_class (),
			   cost_func, cost_data,
			   coarsen_func, coarsen_data,
			   stop_func, stop_data, 
			   fold);
  else
    gts_surface_coarsen (s, 
			 cost_func, cost_data, 
			 coarsen_func, coarsen_data, 
			 stop_func, stop_data, fold);

  /* if verbose on print stats */
  if (verbose) {
    fputc ('\n', stderr);
    gts_surface_print_stats (s, stderr);
    fprintf (stderr, "# volume: %g area: %g\n", 
	     gts_surface_volume (s), gts_surface_area (s));
  }
}

