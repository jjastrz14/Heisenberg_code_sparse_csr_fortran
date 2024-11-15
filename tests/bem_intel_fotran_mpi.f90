! in terminal
!module load intel
!module load impi

program hello_world_mpi 
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror

!This function initializes the MPI environment. It takes in the an error handling variable.
call MPI_INIT(ierror)

!This function returns the total size of the environment in terms of the quantity of processes. 
!The function takes in the MPI environment, an integer to hold the commsize, and an error handling variable.
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)

!This function returns the process id of the process that called the function. 
!The function takes in the MPI environment, an integer to hold the comm rank, and an error handling variable
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

print *, 'Hello World from process: ', process_Rank, 'of ', size_Of_Cluster

!This function cleans up the MPI environment and ends MPI communications.
call MPI_FINALIZE(ierror)

!use mpiifort hello_world_mpi.f90 -o hello_world_mpi.exe

end program hello_world_mpi