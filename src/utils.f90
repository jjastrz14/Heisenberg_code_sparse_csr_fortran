module math_functions
    implicit none
    contains

    subroutine index_map(Sz_subspace_size, ind_1D_target, ind_i, ind_j) 
        implicit none 
        integer(8), intent(in) :: Sz_subspace_size, ind_1D_target
        integer(8), intent(out) :: ind_i, ind_j
        integer(8) :: ind_3
        integer(8) :: ind_Sz_1, ind_Sz_2, ind_tmp, loop_occ, N

        ind_tmp = ind_1D_target
        N = Sz_subspace_size
        loop_occ = 0 

        do while (ind_tmp .GT. 0)
            ind_tmp = ind_tmp - N
            N = N - 1
            loop_occ = loop_occ + 1
        end do 
        ind_i = loop_occ
        ind_j = ind_i + ind_tmp + N 

    end subroutine index_map
    
    subroutine binomialCoefficient(n, k, C)
        integer, intent(in)  :: n, k
        integer(8), intent(out) :: C
        !use, intrinsic :: iso_fortran_env, only: real64, int64
        !integer(int64) :: x
        integer :: i
        integer :: C_temp

        C_temp = 1

        if (k < 0 .or. k > n) then
        C_temp = 0
        else
        do i = 1, min(k, n-k)
            C_temp = C_temp * (n - i + 1) / i
            if (C_temp.GT.2147483647) then
                write(*,*) 'possible integer 32bit 2^32 = 2147483647 overflow'
                error stop
            end if
        end do
        end if

        C = C_temp
    end subroutine binomialCoefficient

    subroutine mmm_mat_mul(D, A, B, C, op_A, op_B, op_C)
        implicit none
        !performs 3 matrix multiplication   D=op(A)*op(B)*op(C)
        double precision, intent(inout) :: D(:,:)
        double precision, intent(in) :: A(:,:), B(:,:), C(:,:)
        CHARACTER*1, intent(in) :: op_A, op_B, op_C
        double precision, allocatable :: temp(:,:)
        integer :: M
        double precision, parameter :: alpha=dcmplx(1.0d0, 0.0d0), beta=dcmplx(0.0d0, 0.0d0)
    
        M=size(D,1)
        allocate( temp(M,M) )
        !matrix multiplication
        !call zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, c, ldc)
        !C := alpha*op(A)*op(B) + beta*C,
        !op(X)=X,X^T,X^H
        !transa, transb - 'n'-op(A)=A, 't'-op(A)=A^T, 'c'-op(A)=A^H-Hermitean conj.
        !alpha and beta are scalars, A, B and C are matrices:
        !A - m-by-k, B - k-by-n, C - m-by-n
        !lda, ldb, ldc - max(1,m)
        call dgemm(op_A, op_B, M, M, M, alpha, A, M, B, M, beta, temp, M)
        call dgemm(op_B, op_C, M, M, M, alpha, temp, M, C, M, beta, D, M)
        deallocate( temp )
    
    end subroutine mmm_mat_mul


    subroutine indexArrayReal(n,Array,Index)
        !credits to https://stackoverflow.com/questions/54860293/ordering-function-in-fortran
        ! this is quick sort from Press Numerical Recipies, subroutine indexx 
        implicit none
        integer, intent(in)  :: n
        double precision   , intent(in)  :: Array(n)
        integer, intent(out) :: Index(n)
        integer, parameter   :: nn=15, nstack=50 ! still don't know why its here
        integer             :: k,i,j,indext,jstack,l,r
        integer             :: istack(nstack)
        double precision             :: a

        if (n.GE.100000000) then 
            write(*,*) 'be carefull for too large input basis sz'
        end if 

        do j = 1,n
            Index(j) = j
        end do
        jstack=0
        l=1
        r=n
        do
            if (r-l < nn) then
                do j=l+1,r
                    indext=Index(j)
                    a=Array(indext)
                    do i=j-1,l,-1
                        if (Array(Index(i)) >= a) exit !<=
                        Index(i+1)=Index(i)
                    end do
                    Index(i+1)=indext
                end do
                if (jstack == 0) return
                r=istack(jstack)
                l=istack(jstack-1)
                jstack=jstack-2
            else
                k=(l+r)/2
                call swap(Index(k),Index(l+1))
                call exchangeIndex(Index(l),Index(r))
                call exchangeIndex(Index(l+1),Index(r))
                call exchangeIndex(Index(l),Index(l+1))
                i=l+1
                j=r
                indext=Index(l+1)
                a=Array(indext)
                do
                    do
                        i=i+1
                        if (Array(Index(i)) <= a) exit
                    end do
                    do
                        j=j-1
                        if (Array(Index(j)) >= a) exit
                    end do
                    if (j > i) exit
                    call swap(Index(i),Index(j))
                end do
                Index(l+1)=Index(j)
                Index(j)=indext
                jstack=jstack+2
                if (jstack > nstack) then
                    write(*,*) 'NSTACK too small in indexArrayReal()'   ! xxx
                    error stop
                end if
                if (r-i+1 <= j-l) then
                    istack(jstack)=r
                    istack(jstack-1)=i
                    r=j-1
                else
                    istack(jstack)=j-1
                    istack(jstack-1)=l
                    l=i
                end if
            end if
        end do
        contains
            subroutine exchangeIndex(i,j)
                integer, intent(inout) :: i,j
                integer              :: swp
                if (Array(j) > Array(i)) then
                    swp=i
                    i=j
                    j=swp
                end if
            end subroutine exchangeIndex
            pure elemental subroutine swap(a,b)
                implicit none
                integer, intent(inout) :: a,b
                integer :: dum
                dum=a
                a=b
                b=dum
            end subroutine swap
    end subroutine indexArrayReal


    subroutine RemoveDuplicates(input_array, final_array)
        integer, dimension(:,:), intent(in) :: input_array
        integer, dimension(:,:), allocatable, intent(out) :: final_array
        integer, dimension(:,:), allocatable :: unique_array
    
        integer :: i, j, k
        logical :: is_duplicate
    
        allocate(unique_array(size(input_array, 1), size(input_array, 2)))
    
        k = 0! Index for the unique_array
    
        do i = 1, size(input_array, 1)
            is_duplicate = .false.
    
            ! Check if the current row is a duplicate of any previous row
            do j = 1, k
                if (all(input_array(i, :) == unique_array(j, :))) then
                    is_duplicate = .true.
                    exit
                end if
            end do
    
            if (.not. is_duplicate) then
                k = k + 1
                unique_array(k, :) = input_array(i, :)
            end if
        end do
         ! Deallocate input_array if it's no longer needed
        ! Resize unique_array to the actual size
        allocate(final_array(k, size(input_array, 2)))
        final_array = unique_array(1:k, :)

    end subroutine RemoveDuplicates


    subroutine FindRowIndex(target_row, search_array, index)
        integer, dimension(:), intent(in) :: target_row
        integer, dimension(:,:), intent(in) :: search_array
        integer, intent(out):: index
        integer :: i
    
            ! Initialize the function value to indicate not found
            index = 0
        
            do i = 1, size(search_array, 1)
                if (all(target_row == search_array(i, :))) then
                    index = i
                    exit ! Exit the loop when found
                end if
            end do

    end subroutine FindRowIndex

end module math_functions


module timing_utilities
    implicit none
    private              ! Make everything private by default
    public :: timer      ! Only expose what we want to use
    
    ! Define the supported time units as enumerated type
    type, public :: time_unit_type
        integer :: seconds = 1
        integer :: minutes = 2
        integer :: hours = 3
    end type time_unit_type
    
    type(time_unit_type), parameter, public :: time_unit = time_unit_type()
    
    ! Constants for time conversion
    double precision, parameter :: SECONDS_TO_MINUTES = 1.0d0/60.0d0
    double precision, parameter :: SECONDS_TO_HOURS = 1.0d0/3600.0d0
    
    ! Timer type definition with all necessary components
    type timer
        private
        integer(8) :: start_count = 0       ! Start timestamp
        integer(8) :: end_count = 0         ! End timestamp
        integer(8) :: count_rate            ! System clock rate
        double precision :: total_time = 0.0d0  ! Accumulated time
        logical :: is_running = .false.     ! Timer state flag
        contains
            ! Timer methods
            procedure :: start => start_timer
            procedure :: stop => stop_timer
            procedure :: reset => reset_timer
            procedure :: get_elapsed => get_elapsed_time
            procedure :: print_elapsed => print_elapsed_time
    end type timer
    
    contains
        ! Start the timer if it's not already running
        subroutine start_timer(this)
            class(timer), intent(inout) :: this
            
            if (.not. this%is_running) then
                call system_clock(count_rate=this%count_rate)
                call system_clock(this%start_count)
                this%is_running = .true.
            end if
        end subroutine start_timer
        
        ! Stop the timer and accumulate elapsed time
        subroutine stop_timer(this)
            class(timer), intent(inout) :: this
            
            if (this%is_running) then
                call system_clock(this%end_count)
                this%total_time = this%total_time + &
                    dble(this%end_count - this%start_count) / dble(this%count_rate)
                this%is_running = .false.
            end if
        end subroutine stop_timer
        
        ! Reset the timer to initial state
        subroutine reset_timer(this)
            class(timer), intent(inout) :: this
            
            this%total_time = 0.0d0
            this%is_running = .false.
        end subroutine reset_timer
        
        ! Get elapsed time in requested units
        function get_elapsed_time(this, unit) result(elapsed)
            class(timer), intent(in) :: this
            integer, intent(in), optional :: unit
            double precision :: elapsed
            
            if (.not. present(unit)) then
                elapsed = this%total_time    ! Default to seconds
                return
            end if
            
            select case(unit)
                case(time_unit%seconds)
                    elapsed = this%total_time
                case(time_unit%minutes)
                    elapsed = this%total_time * SECONDS_TO_MINUTES
                case(time_unit%hours)
                    elapsed = this%total_time * SECONDS_TO_HOURS
                case default
                    elapsed = this%total_time
            end select
        end function get_elapsed_time
        
        ! Print elapsed time with optional unit and label
        subroutine print_elapsed_time(this, unit, label)
            class(timer), intent(in) :: this
            integer, intent(in), optional :: unit
            character(len=*), intent(in), optional :: label
            double precision :: elapsed_time
            character(len=:), allocatable :: time_unit_str, prefix
            
            ! Set the label (default or custom)
            if (present(label)) then
                prefix = label
            else
                prefix = "Elapsed time"
            endif
            
            ! Get the time in requested units
            elapsed_time = this%get_elapsed(unit)
            
            ! Determine the appropriate unit string
            if (present(unit)) then
                select case(unit)
                    case(time_unit%seconds)
                        time_unit_str = "seconds"
                    case(time_unit%minutes)
                        time_unit_str = "minutes"
                    case(time_unit%hours)
                        time_unit_str = "hours"
                    case default
                        time_unit_str = "seconds"
                end select
            else
                time_unit_str = "seconds"
            endif
            
            ! Print with consistent formatting
            write(*,'(1X,A,": ",F15.6," ",A)') prefix, elapsed_time, time_unit_str
            
        end subroutine print_elapsed_time
end module timing_utilities