#include "TA_common.h"

double best_of_rounded_delta(int *xn_plus)
{
  double fxc;
  int j, d = n_dimensions;
  int xn_plus_grow[d];

  // Growing, shrinking.
  // Grower, shrinker that copy the point
  for (j = 0; j < d; j++)
    xn_plus_grow[j] = xn_plus[j];
  grow_box_randomly(xn_plus_grow);

  // Now, create the official numbers.
  // official update from modified points
  fxc = get_delta(xn_plus_grow);

  return fxc;
}

void update_points(double* fxc, double* current, int* xc_index, int* xn_pm_index, int* xn_extraminus_index, int* k, int mc) {
      // generation of random point xc
      int* xn_plus_index = xn_pm_index;
      generate_xc_delta(xc_index);

      //(Possibly) Snaps the point upwards and computes the fitness
      *current = best_of_rounded_delta(xc_index);

      // draw a neighbour of xc
      generate_neighbor_delta(xn_plus_index, xc_index, k, mc);

      // Compute the threshold
      *fxc = best_of_rounded_delta(xn_plus_index);
}

void generate_neighbor_deltabardelta(int* xc_index, int* xn_pm_index, int* xn_extraminus_index, int* k, int mc) {
  generate_neighbor_delta(xn_pm_index, xc_index, k, mc);
};

double best_of_rounded_deltabardelta(int* xc_index, int* xn_pm_index, int* xn_extraminus_index, int* xn_best_index) {
  return best_of_rounded_delta(xn_pm_index);
}

void generate_xc_deltabardelta(int* xc_index, int* xn_pm_index, int* xn_extraminus_index) {
  generate_xc_delta(xc_index);
}
