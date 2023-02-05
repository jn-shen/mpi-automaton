# MPP Coursework Assessment
## File Structure
* automaton.c : main function
* mpi_utils.c : mpi-related functions
    * cartesian_create: 
    * split_grid
    * init_halos_cell
    * scatter_cell
    * gather_cell
    * prepare_swap_halos
    * swap_halos
* game_utils.c : game-related functions
    * init_grid
    * update_cell
    * terminate_condition
* cellio.c : visualisation
    * cellwrite
    * cellwritedynamic
* unirand.c : random numbers
    * uni
    * rstart
    * rinit
* automaton.h : header files
* automaton.job : slurm job script
* Makefile 
## Compiling
You must load these modules in order to run this program
```
module load mpt
module load intel-compilers-19 
```
You can use  ```make``` to compile this program
## Running

### Command Line

```
mpirun -n 16 ./automaton 1234
``` 
### Batch System
```
sbatch automaton.job
```
## display image
```
module load ImageMagick
display cell.pbm
```

## Notes
* You can change the size of L in automaton.h
* 16 is the number of processes, you can change this value to whatever you want except 0 or the number of processes larger than L*L.  
In automaton.job, you can see:
```
srun -n 16 ./automaton 1234
```
* if you want to change the number of processes in automaton.job, for every 36 processes you need to add a node. You can change this ```#SBATCH --nodes``` to add.
* 1234 is a random seed, which is necessary for the program to run. Users can change this value.
