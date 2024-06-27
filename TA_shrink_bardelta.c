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

double oldmain(struct grid *grid, double **pointset, int n, int d, int mc, int i_tilde, int trials)
{
  int k[d], start[d];

  int i, j, p, t; // loop variables

  double thresh[i_tilde + 1]; // Thresholdsequence
  double T;               // current Threshold

  double fxc;
  int xc_index[d], xn_minus_index[d], xn_extraminus_index[d];
  int xn_best_index[d]; // Indices of current point, neighbour
  double current, global[trials + 1], best; // current and global best values

  int outerloop = i_tilde, innerloop = i_tilde;

  // Sort the grid points, setup global variables
  process_coord_data(grid, pointset, n, d);

  // Algorithm starts here
  for (t = 1; t <= trials; t++)
  { // Initialization
    // Initialize k-value
    for (j = 0; j < d; j++)
    {
      start[j] = (int)((grid->n_coords[j] - 1) / 2);
    }
    // Initialize mc-value
    mc = 2;

    // Initialize iteration count
    int current_iteration = 0;

    // Generate threshold sequence
    for (i = 1; i <= outerloop; i++)
    {
      current_iteration++;
      // Update k-value
      for (j = 0; j < d; j++)
      {
        k[j] = start[j] * (((double)outerloop - current_iteration) / (outerloop)) +
               1 * ((double)current_iteration / (outerloop));
        //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
      }

      // Update mc-value
      mc = 2 + (int)(current_iteration / outerloop * (d - 2));

      // generation of random point xc
      generate_xc_bardelta(grid, xn_minus_index, xn_extraminus_index);

      //(Possibly) Snap the points and compute the largest of the rounded values
      current = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);

      // draw a neighbour of xc
      generate_neighbor_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index, k, mc);

      // Compute the threshold
      fxc = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index);
      thresh[i] = 0.0 - fabs(fxc - current);
    }

    // sort the thresholds in increasing order
    quicksort(1, outerloop, thresh);

    current = 0;
    global[t] = 0;

    // Initialize k-value
    for (j = 0; j < d; j++)
    {
      start[j] = (int)((grid->n_coords[j] - 1) / 2);
    }

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
        for (j = 0; j < d; j++)
        {
          k[j] = start[j] * (((double)innerloop * outerloop - current_iteration) / (innerloop * outerloop)) +
                 1 * ((double)current_iteration / (innerloop * outerloop));
          //		k[j]=(int)(start[j]-(int)(current_iteration/(innerloop*outerloop)*(start[j]-1)));
        }

        // Update mc-value
        mc = 2 + (int)(current_iteration / (innerloop * outerloop) * (d - 2));
        // mc=2;

        // Get random neighbor
        generate_neighbor_bardelta(grid, xn_minus_index, xn_extraminus_index, xc_index, k, mc);

        //(Possibly) Snap the points and compute the best of the rounded points
        fxc = best_of_rounded_bardelta(grid, xn_minus_index, xn_extraminus_index, xn_best_index);
        // Global update if necessary
        if (fxc > global[t])
        {
          global[t] = fxc;
        }
        // Update of current best value if necessary
        if (fxc - current >= T)
        {
          current = fxc;
          for (j = 0; j < d; j++)
          {
            xc_index[j] = xn_best_index[j];
          }
        }
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
  printf("%g\n", oldmain(&grid, param.pointset, param.npoints, param.dim, param.mc, param.i_tilde, param.trials));
  return EXIT_SUCCESS;
}
