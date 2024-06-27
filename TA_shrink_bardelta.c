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
  int i, j, p, t; // loop variables

  double thresh[i_tilde + 1]; // Thresholdsequence
  double T;               // current Threshold

  int xc_index[d], xn_minus_index[d], xn_extraminus_index[d];
  int xn_best_index[d]; // Indices of current point, neighbour
  double current, global[trials + 1], best; // current and global best values
  struct kmc kmc;
  int _k[d];
  kmc.k = _k;
  int outerloop = i_tilde, innerloop = i_tilde;

  // Sort the grid points, setup global variables
  process_coord_data(grid, pointset, n, d);

  // Algorithm starts here
  for (t = 1; t <= trials; t++)
  { // Initialization
    // Initialize iteration count
    int current_iteration = 0;

    // Generate threshold sequence
    for (i = 1; i <= outerloop; i++)
    {
      current_iteration++;
      // Update k-value
      // Update mc-value
      get_kmc(grid, &kmc, current_iteration, outerloop);
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
    quicksort(1, outerloop, thresh);

    current = 0;
    global[t] = 0;
    // draw a random initial point
    generate_xc_bardelta(grid, xn_minus_index, xn_extraminus_index);

    //(Possibly) Snap and compute the best of the rounded points and update current value
    current = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);

    global[t] = current;

    current_iteration = 0;
    for (i = 1; i <= outerloop; i++)
    {
      T = thresh[i];

      for (p = 1; p <= innerloop; p++)
      {
        current_iteration++;

        // Update k-value
        // Update mc-value
        get_kmc(grid, &kmc, current_iteration, innerloop * outerloop);

        // Get random neighbor
        generate_neighbor_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index, kmc.k, kmc.mc);

        //(Possibly) Snap the points and compute the best of the rounded points
        double fxc = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xn_best_index);
        if (fxc > global[t])
        {
          global[t] = fxc;
        }
        ta_update_point(fxc, &current, T, xc_index, xn_best_index, d);
      } // innerloop
    } // outerloop
  } // trials

  // best calculated value
  best = 0;
  for (t = 1; t <= trials; t++)
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
  return EXIT_SUCCESS;
}
