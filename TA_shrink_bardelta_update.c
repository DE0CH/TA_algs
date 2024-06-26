#include "TA_common.h"

double best_of_rounded_bardelta(int *xn_minus, int *xn_extraminus, int *xc_index)
{
  double fxn_extraminus;
  double fxc;
  int j, d = n_dimensions;
  int use_extraminus = 0;
  int xn_minus_snap[d], xn_extraminus_snap[d];
  for (j = 0; j < d; j++)
    if (xn_minus[j] != xn_extraminus[j])
    {
      use_extraminus = 1;
      break;
    }

  // Growing, shrinking.
  // Grower, shrinker that copy the point
  for (j = 0; j < d; j++)
    xn_minus_snap[j] = xn_minus[j];
  snap_box(xn_minus_snap);
  if (use_extraminus)
  {
    for (j = 0; j < d; j++)
      xn_extraminus_snap[j] = xn_extraminus[j];
    snap_box(xn_extraminus_snap);
  }

  // Now, create the official numbers.
  // official update from modified points
  fxc = get_bar_delta(xn_minus_snap);
  if (use_extraminus)
  {
    fxn_extraminus = get_bar_delta(xn_extraminus_snap);
    fxc = max(fxc, fxn_extraminus);
  }

  // Remains only to copy the winning point to output variable xc_index.
  if (use_extraminus && (fxn_extraminus >= fxc))
  {
    for (j = 0; j < d; j++)
      xc_index[j] = xn_extraminus[j];
  }
  else
  {
    for (j = 0; j < d; j++)
      xc_index[j] = xn_minus[j];
  }

  return fxc;
}

void update_points(double* fxc, double* current, int* xc_index, int* xn_pm_index, int* xn_extraminus_index, int* k, int mc) {
  // generation of random point xc
  int* xn_minus_index = xn_pm_index;
  generate_xc_bardelta(xn_minus_index, xn_extraminus_index);

  //(Possibly) Snap the points and compute the largest of the rounded values
  *current = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);

  // draw a neighbour of xc
  generate_neighbor_bardelta(xn_minus_index, xn_extraminus_index, xc_index, k, mc);
  // Compute the threshold
  *fxc = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);
}

void generate_neighbor_deltabardelta(int* xc_index, int* xn_pm_index, int* xn_extraminus_index, int* k, int mc) {
  generate_neighbor_bardelta(xn_pm_index, xn_extraminus_index, xc_index, k, mc);
};

double best_of_rounded_deltabardelta(int* xn_pm_index, int* xn_extraminus_index, int* xn_best_index) {
  return best_of_rounded_bardelta(xn_pm_index, xn_extraminus_index, xn_best_index);
}


void generate_xc_deltabardelta(int* xn_pm_index, int* xn_extraminus_index) {
  generate_xc_bardelta(xn_pm_index, xn_extraminus_index);
}
