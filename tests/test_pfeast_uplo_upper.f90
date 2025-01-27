program hello_world_pfeast_local
    implicit none
    include 'mpif.h'
    ! 4x4 global eigenvalue system == two 2x4 local matrices
    integer, parameter :: Nloc=2, NNZloc=7
    character(len=1) :: UPLO='U'
    double precision, allocatable :: Aloc(:)
    integer, allocatable :: IAloc(:)
    integer, allocatable :: JAloc(:)
    ! Input parameters for FEAST
    integer, dimension(64) :: fpm
    integer :: M0=3 ! Search subspace dimension
    double precision :: Emin=3.0d0, Emax=5.0d0 ! Search interval
    ! Output variables for FEAST
    double precision, dimension(:), allocatable :: E, res
    double precision, dimension(:,:), allocatable :: X
    double precision :: epsout
    integer :: nL3, rank3, loop, info, M, i
    ! MPI
    integer :: code

    call MPI_INIT(code)

    ! Allocate memory for eigenvalues, eigenvectors, residual
    allocate(E(M0), X(Nloc,M0), res(M0))

    ! Initialize PFEAST and distribute matrix
    nL3 = 2
    call pfeastinit(fpm, MPI_COMM_WORLD, nL3)

    call MPI_COMM_RANK(fpm(49), rank3, code) ! Find rank of new L3 communicator

    if (rank3 == 0) then
        allocate(Aloc(6), IAloc(3), JAloc(6))
        Aloc = (/2.0d0, -1.0d0, -1.0d0, 3.0d0, -1.0d0, -1.0d0/)
        IAloc = (/1, 4, 7/)
        JAloc = (/1, 2, 3, 2, 3, 4/)
    else if (rank3 == 1) then
        allocate(Aloc(3), IAloc(3), JAloc(3))
        Aloc = (/3.0d0, -1.0d0, 2.0d0/)
        IAloc = (/1, 3, 4/)
        JAloc = (/3, 4, 4/)
    endif

    ! PFEAST
    fpm(1) = 1 ! Change from default value (print info on screen)
    call pdfeast_scsrev(UPLO, Nloc, Aloc, IAloc, JAloc, fpm, epsout, loop, Emin, Emax, M0, E, X, M, res, info)

    ! Report results
    if (info == 0) then
        print *, 'Solutions (Eigenvalues/Eigenvectors/Residuals) at rank L3', rank3
        do i = 1, M
            print *, 'Eigenvalue', i
            print *, 'E =', E(i), 'X =', X(:, i), 'Res =', res(i)
        enddo
    endif

    deallocate(E, X, res)
    deallocate(Aloc, IAloc, JAloc)
    call MPI_FINALIZE(code)
end program hello_world_pfeast_local