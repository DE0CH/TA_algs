/*
   reads a point set and computes a lower bound on its star discrepancy
*/

/*
  i_tilde inner and outer loops (effect of alpha is considered to be negligible)
  threshold aus Paaren von Nachbarn berechnet
*/

// This time, the file is written so that first element is x[0] :-p

// compute threshold sequence once, then reuse it
// not using alpha

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "TA_common.h"

double oldmain(struct grid *grid, double **pointset, int n, int d, int i_tilde, int trials)
{
  int search_population = 200;
  double thresh[i_tilde]; // Thresholdsequence

  int xc_index[d], xn_minus_index[d], xn_extraminus_index[d];
  int xn_best_index[d]; // Indices of current point, neighbour
  double current, global[trials], best; // current and global best values
  struct kmc kmc;
  int _k[d];
  kmc.k = _k;
  struct history history = init_history(d, trials*i_tilde*i_tilde + trials);
  double search_points[search_population * d];
  populate_random_search_points(search_population, d, search_points);

  // Sort the grid points, setup global variables
  process_coord_data(grid, pointset, n, d);

  // Algorithm starts here
  for (int t = 0; t < trials; t++)
  { // Initialization
    // Initialize iteration count

    // Generate threshold sequence
    for (int i = 0; i < i_tilde; i++)
    {
      int current_iteration = i + 1;
      // Update k-value
      // Update mc-value
      get_kmc(grid, &kmc, current_iteration, i_tilde);
      // generation of random point xc
      generate_xc_bardelta(grid, xn_minus_index, xn_extraminus_index);

      //(Possibly) Snap the points and compute the largest of the rounded values
      current = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);

      // draw a neighbour of xc
      generate_neighbor_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index, kmc.k, kmc.mc);

      // Compute the threshold
      double fxc = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);
      thresh[i] = 0.0 - fabs(fxc - current);
    }

    // sort the thresholds in increasing order
    quicksort(0, i_tilde, thresh);

    current = 0;
    global[t] = 0;
    // draw a random initial point
    if (t == 0) {
      generate_xc_bardelta(grid, xn_minus_index, xn_extraminus_index);
    } else {
      // populate_random_search_points(search_population, d, search_points);
      double pp[d];
      furthest_point(&history, search_population, search_points, pp);
      record_history_raw(&history, pp);
      round_point_down(grid, pp, xn_minus_index);
      round_point_extradown(grid, pp, xn_extraminus_index);
    }

    current = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);

    //(Possibly) Snap and compute the best of the rounded points and update current value

    global[t] = current;

    for (int i = 0; i < i_tilde; i++)
    {
      double T = thresh[i];

      for (int p = 0; p < i_tilde; p++)
      {
        int current_iteration = i*i_tilde + p + 1;

        // Update k-value
        // Update mc-value
        get_kmc(grid, &kmc, current_iteration, i_tilde * i_tilde);

        // Get random neighbor
        generate_neighbor_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index, kmc.k, kmc.mc);

        //(Possibly) Snap the points and compute the best of the rounded points
        double fxc = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xn_best_index);
        if (fxc > global[t])
        {
          global[t] = fxc;
        }
        ta_update_point(fxc, &current, T, xc_index, xn_best_index, d);
        record_history(grid, &history, xn_best_index);
      } // innerloop
    } // outerloop
  } // trials

  free_history(&history);
  // best calculated value
  best = 0;
  for (int t = 0; t < trials; t++)
  {
    if (global[t] > best)
    {
      best = global[t];
    }
  }
  return best;
}

int main(int argc, char **argv)
{
  struct grid grid;
  struct initial_params param;
  fprintf(stderr, "Calling Carola calculation\n");
  read_points(argc, argv, &param);
  printf("%g\n", oldmain(&grid, param.pointset, param.npoints, param.dim, param.i_tilde, param.trials));
  free_initial_params(&param);
  free_grid(&grid);
  return EXIT_SUCCESS;
}
