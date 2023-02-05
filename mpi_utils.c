#include "automaton.h"

int stride_x, stride_y; // standard stride used for scatter and gather to locate

/***********   cartesian_create   ***********/
/*
 *  Functionality: Create a Cartesian topology
 *  Parameters:
 *  old_comm: Input communicator
 *  size: Number of processes
 *  dims: Number of processes in each dimension
 *  new_comm: Communicator with new cartesian topology
 */
void cartesian_create(MPI_Comm old_comm, int size, int dims[2], MPI_Comm *new_comm)
{
  int rank;
  int reorder = FALSE;
  // The second dimension is periodic
  int period[2] = {FALSE, TRUE};
  dims[0] = 0;
  dims[1] = 0;
  // Generate dimension
  MPI_Dims_create(size, 2, dims);
  // Create cartesian topology
  MPI_Cart_create(old_comm, 2, dims, period, reorder, new_comm);
}

/***********   split_grid   ***********/
/*
 *  Functionality: 2D decomposite the grid of cells for each process, it can also handle not decomposing exactly
 *  Parameters:
 *  rank: Rank of the calling process
 *  dims: Number of processes in each dimension
 *  comm: Input communicator
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 */
void split_grid(int rank, int dims[2], MPI_Comm comm, int *LX, int *LY)
{
  int coords[2];
  int remainder_x = L % dims[0];
  int remainder_y = L % dims[1];
  *LX = L / dims[0];
  *LY = L / dims[1];
  stride_x = *LX;
  stride_y = *LY;
  // get the cartesian coordinate
  MPI_Cart_coords(comm, rank, 2, coords);
  // The case that not decomposing exactly
  if (coords[0] == dims[0] - 1)
  {
    *LX = *LX + remainder_x;
  }
  if (coords[1] == dims[1] - 1)
  {
    *LY = *LY + remainder_y;
  }
}

/***********   init_halos_cell   ***********/
/*
 *  Functionality: Add halos to the cells assigned to each process
 *  Parameters:
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  small_cell: Cells without halos assigned to each process
 *  cell: Cells with halos assigned to each process
 */
void init_halos_cell(int LX, int LY, int small_cell[LX][LY], int cell[LX + 2][LY + 2])
{
  int i;
  int j;
  for (i = 1; i <= LX; i++)
  {
    for (j = 1; j <= LY; j++)
    {
      cell[i][j] = small_cell[i - 1][j - 1];
    }
  }
  for (i = 0; i <= LX + 1; i++)
  {
    cell[i][0] = 0;
    cell[i][LY + 1] = 0;
  }

  for (j = 0; j <= LY + 1; j++)
  {
    cell[0][j] = 0;
    cell[LX + 1][j] = 0;
  }
}

/***********   scatter_cell   ***********/
/*
 *  Functionality: Assign cells to each process
 *  Parameters:
 *  rank: Rank of the calling process
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  all_cell: All cells for the grid
 *  comm: Input communicator
 *  small_cell: Cells without halos assigned to each process
 */
void scatter_cell(int rank, int LX, int LY, int all_cell[L][L], MPI_Comm comm, int small_cell[LX][LY])
{
  int coords[2];
  // broadcast to all processes
  MPI_Bcast(&all_cell[0][0], L * L, MPI_INT, 0, comm);
  // get the cartesian coordinate
  MPI_Cart_coords(comm, rank, 2, coords);
  // get the parts needed by each process
  for (int i = 0; i < LX; i++)
  {
    for (int j = 0; j < LY; j++)
    {
      small_cell[i][j] = all_cell[coords[0] * stride_x + i][coords[1] * stride_y + j];
    }
  }
}

/***********   gather_cell   ***********/
/*
 *  Functionality: Gather cells from each process
 *  Parameters:
 *  rank: Rank of the calling process
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  small_cell: Cells without halos assigned to each process
 *  cell: Cells with halos assigned to each process
 *  comm: Input communicator
 *  all_cell: All cells for the grid
 */
void gather_cell(int rank, int LX, int LY, int small_cell[LX][LY], int cell[LX + 2][LY + 2], MPI_Comm comm, int all_cell[L][L])
{
  int i, j;
  int coords[2];
  int tmpcell[L][L];
  // remove halos
  for (i = 1; i <= LX; i++)
  {
    for (j = 1; j <= LY; j++)
    {
      small_cell[i - 1][j - 1] = cell[i][j];
    }
  }
  // initialize the temporary array
  for (i = 0; i < L; i++)
  {
    for (j = 0; j < L; j++)
    {
      tmpcell[i][j] = 0;
    }
  }
  // assign values to the temporary array in the correct position
  MPI_Cart_coords(comm, rank, 2, coords);
  for (i = 0; i < LX; i++)
  {
    for (j = 0; j < LY; j++)
    {
      tmpcell[coords[0] * stride_x + i][coords[1] * stride_y + j] = small_cell[i][j];
    }
  }
  MPI_Reduce(&tmpcell[0][0], &all_cell[0][0], L * L, MPI_INT, MPI_SUM, 0, comm);
}

/***********   prepare_swap_halos   ***********/
/*
 *  Functionality: Define the important parameter for halos swapping
 *  Parameters:
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  comm: Input communicator
 *  left, right, top, down: The four processes that adjoin this process
 *  vector: Defined vector datatype
 */
void prepare_swap_halos(int LX, int LY, MPI_Comm comm, int *left, int *right, int *top, int *down, MPI_Datatype *vector)
{
  int disp = 1;
  MPI_Type_vector(LX, 1, LY + 2, MPI_INT, vector);
  MPI_Type_commit(vector);
  MPI_Cart_shift(comm, 0, disp, left, right);
  MPI_Cart_shift(comm, 1, disp, down, top);
}

/***********   swap_halos   ***********/
/*
 *  Functionality: Swap halos between processes
 *  Parameters:
 *  left, right, top, down: The four processes that adjoin this process
 *  LX: The local dimension of the cell in the x-axis direction for each process
 *  LY: The local dimension of the cell in the y-axis direction for each process
 *  vector: Defined vector datatype
 *  comm: Input communicator
 *  cell: Cells with halos assigned to each process
 */
void swap_halos(int left, int right, int top, int down, int LX, int LY, MPI_Datatype vector, MPI_Comm comm, int cell[LX + 2][LY + 2])
{
  MPI_Request request[4];
  MPI_Status status[4];
  // vertical halos send 
  MPI_Issend(&cell[1][1], LY, MPI_INT, left, 0, comm, &request[0]);
  MPI_Issend(&cell[LX][1], LY, MPI_INT, right, 0, comm, &request[1]);
  // horizontal halos send
  MPI_Issend(&cell[1][1], 1, vector, down, 0, comm, &request[2]);
  MPI_Issend(&cell[1][LY], 1, vector, top, 0, comm, &request[3]);
  // receive halos
  MPI_Recv(&cell[LX + 1][1], LY, MPI_INT, right, 0, comm, &status[0]);
  MPI_Recv(&cell[0][1], LY, MPI_INT, left, 0, comm, &status[1]);
  MPI_Recv(&cell[1][LY + 1], 1, vector, top, 0, comm, &status[2]);
  MPI_Recv(&cell[1][0], 1, vector, down, 0, comm, &status[3]);
  // wait all processes finish
  MPI_Waitall(4, request, status);
}