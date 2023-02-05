#include <mpi.h>

//System size L
#define L 768
#define TRUE 1
#define FALSE 0


//Initialize the grid of cells using the specified random number seed
int init_grid(char *str_seed[], int size, int all_cell[L][L], int max_step);
//Update the state of the cells and calculate the number of local live cells and the number of local cells that changed their values
int update_cell(int LX, int LY, MPI_Comm comm, int *local_changed_cell, int cell[LX + 2][LY + 2]);
//Check if the termination condition is met and calculate the number of live cells and the number of cells that changed their values
int terminate_condition(int num_local_cell, int local_changed_cell, int *num_cell, int *changed_cell, int max_num_cell, int min_num_cell, MPI_Comm comm);

//Create a Cartesian topology
void cartesian_create(MPI_Comm old_comm, int size, int dims[2], MPI_Comm *new_comm);
//2D decomposite the grid of cells for each process
void split_grid(int rank, int dims[2], MPI_Comm comm, int *LX, int *LY);
//Add halos to the cells assigned to each process
void init_halos_cell(int LX, int LY, int small_cell[LX][LY], int cell[LX + 2][LY + 2]);
//Assign cells to each process
void scatter_cell(int rank, int LX, int LY, int all_cell[L][L], MPI_Comm comm, int small_cell[LX][LY]);
//Gather cells from each process
void gather_cell(int rank, int LX, int LY, int small_cell[LX][LY], int cell[LX + 2][LY + 2], MPI_Comm comm, int all_cell[L][L]);
//Define the important parameter for halos swapping
void prepare_swap_halos(int LX, int LY, MPI_Comm comm, int *left, int *right, int *top, int *down, MPI_Datatype *vector);
//Swap halos between processes
void swap_halos(int left, int right, int top, int down, int LX, int LY, MPI_Datatype vector, MPI_Comm comm, int cell[LX + 2][LY + 2]);

//Visualisation
void cellwrite(char *cellfile, int cell[L][L]);
void cellwritedynamic(char *cellfile, int **cell, int l);

//Random numbers
void rinit(int ijkl);
float uni(void);
