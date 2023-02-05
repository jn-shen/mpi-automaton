#include <stdio.h>
#include <stdlib.h>
#include "automaton.h"

/***********   init_grid   ***********/
/*
 *  Functionality: Initialize the grid of cells using the specified random number seed
 *  Parameters:
 *  str_seed: Random seed for char type
 *  size: Number of processes
 *  all_cell: All cells for the map
 *  max_step: The maximum execution steps of the game
 *
 *  Output:
 *  num_cell: Number of live cells
 */
int init_grid(char *str_seed[], int size, int all_cell[L][L], int max_step)
{
    int seed;
    double rho;
    double r;
    int num_cell;
    printf("automaton: running on %d process(es)\n", size);
    rho = 0.49;
    seed = atoi(str_seed[1]);
    printf("automaton: L = %d, rho = %f, seed = %d, max_step = %d\n", L, rho, seed, max_step);
    rinit(seed);
    num_cell = 0;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            r = uni();
            if (r < rho)
            {
                all_cell[i][j] = 1;
                num_cell = num_cell + 1;
            }
            else
            {
                all_cell[i][j] = 0;
            }
        }
    }
    printf("automaton: rho = %f, live cells = %d, actual density = %f\n", rho, num_cell, ((double)num_cell) / ((double)L * L));
    return num_cell;
}

/***********   update_cell   ***********/
/*
 *  Functionality: Update the state of the cells and calculate the number of local live cells and the number of local cells that changed their values
 *  Parameters:
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  comm: Input communicator
 *  local_changed_cell: Number of changed cells for each process
 *  cell: Cells with halos calculated by each process
 *
 *  Output:
 *  num_local_cell: Number of live cells for each process
 */
int update_cell(int LX, int LY, MPI_Comm comm, int *local_changed_cell, int cell[LX + 2][LY + 2])
{
    int i, j, num_local_cell;
    int neighbour[LX + 2][LY + 2];
    int changed_local_cell;
    for (i = 1; i <= LX; i++)
    {
        for (j = 1; j <= LY; j++)
        {
            neighbour[i][j] = cell[i][j] + cell[i][j - 1] + cell[i][j + 1] + cell[i - 1][j] + cell[i + 1][j];
        }
    }
    num_local_cell = 0;
    changed_local_cell = 0;
    for (i = 1; i <= LX; i++)
    {
        for (j = 1; j <= LY; j++)
        {
            // game rules
            if (neighbour[i][j] == 2 || neighbour[i][j] == 4 || neighbour[i][j] == 5)
            {
                // if cell is dead
                if (cell[i][j] == 0)
                {
                    changed_local_cell++;
                }
                cell[i][j] = 1;
                num_local_cell++;
            }
            else
            {
                // if cell is live
                if (cell[i][j] == 1)
                {
                    changed_local_cell++;
                }
                cell[i][j] = 0;
            }
        }
    }
    *local_changed_cell = changed_local_cell;
    return num_local_cell;
}

/***********   terminate_condition   ***********/
/*
 *  Functionality: Check if the termination condition is met and calculate the number of live cells and the number of cells that changed their values
 *  Parameters:
 *  num_local_cell: Number of live cells for each process
 *  local_changed_cell: Number of changed cells for each process
 *  num_cell: Number of live cells
 *  changed_cell: Number of changed cells
 *  max_num_cell: Maximum number of viable cells
 *  min_num_cell: Minimum number of viable cells
 *  comm: Input communicator
 *
 *  Output:
 *  1 represents the termination condition is met
 *  0 represents the not met termination condition
 */
int terminate_condition(int num_local_cell, int local_changed_cell, int *num_cell, int *changed_cell, int max_num_cell, int min_num_cell, MPI_Comm comm)
{
    // reduce to rank 0
    MPI_Reduce(&local_changed_cell, changed_cell, 1, MPI_INT, MPI_SUM, 0, comm);
    // reduce to all processes
    MPI_Allreduce(&num_local_cell, num_cell, 1, MPI_INT, MPI_SUM, comm);
    // terminate conditon
    if (*num_cell > max_num_cell || *num_cell < min_num_cell)
    {
        return 1;
    }
    else
        return 0;
}