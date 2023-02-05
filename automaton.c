#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

int main(int argc, char *argv[])
{
  int LX, LY, size, rank, num_cell, max_num_cell, min_num_cell, max_step, print_freq, left, right, top, down, num_local_cell, local_changed_cell, changed_cell;
  int dims[2], all_cell[L][L];
  double time_start, time_stop;
  MPI_Comm comm, comm_2d;
  MPI_Datatype horizontal_halos;
  // Initialize the MPI program and get rank and size
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  // Create cartesian topology
  cartesian_create(comm, size, dims, &comm_2d);
  // 2D decomposite the grid of cells
  split_grid(rank, dims, comm_2d, &LX, &LY);
  // Define array without halos and with halos
  int small_cell[LX][LY], cell[LX + 2][LY + 2];
  // Initialize the grid of cells and get the number of live cells
  max_step = 10 * L;
  print_freq = 500;
  if (rank == 0)
  {
    num_cell = init_grid(argv, size, all_cell, max_step);
  }
  // Set terminate conditions for each process
  MPI_Bcast(&num_cell, 1, MPI_INT, 0, comm);
  max_num_cell = num_cell * 3 / 2;
  min_num_cell = num_cell * 2 / 3;
  // Scatter cells to each process
  scatter_cell(rank, LX, LY, all_cell, comm_2d, small_cell);
  // Add halos to cells
  init_halos_cell(LX, LY, small_cell, cell);
  // Define the important paramater for halos swapping
  prepare_swap_halos(LX, LY, comm_2d, &left, &right, &top, &down, &horizontal_halos);
  // Start the timer
  MPI_Barrier(MPI_COMM_WORLD);
  time_start = MPI_Wtime();
  // Start the game
  for (int step = 1; step <= max_step; step++)
  {
    swap_halos(left, right, top, down, LX, LY, horizontal_halos, comm_2d, cell);
    num_local_cell = update_cell(LX, LY, comm_2d, &local_changed_cell, cell);
    if (terminate_condition(num_local_cell, local_changed_cell, &num_cell, &changed_cell, max_num_cell, min_num_cell, comm_2d))
    {
      break;
    }
    if (step % print_freq == 0)
    {
      if (rank == 0)
      {
        printf("automaton: number of live cells on step %d is %d, ", step, num_cell);
        printf("%d cells change their value\n", changed_cell);
      }
    }
  }
  // Stop the timer
  MPI_Barrier(MPI_COMM_WORLD);
  time_stop = MPI_Wtime();
  // Gather cells to each process
  gather_cell(rank, LX, LY, small_cell, cell, comm_2d, all_cell);
  // Write image and print execute time
  if (rank == 0)
  {
    cellwrite("cell.pbm", all_cell);
    printf("time:%f\n", time_stop - time_start);
  }
  MPI_Finalize();
  return 0;
}
