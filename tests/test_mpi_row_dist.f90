subroutine distribute_rows(N, num_processes, start_rows, end_rows, elements, ierr)
  implicit none
  integer, intent(in) :: N                ! Size of the matrix (N x N)
  integer, intent(in) :: num_processes    ! Number of MPI processes
  integer, intent(out) :: start_rows(num_processes) ! Start row for each process
  integer, intent(out) :: end_rows(num_processes)   ! End row for each process
  integer, intent(out) :: elements(num_processes)   ! Number of elements per process
  integer, intent(out) :: ierr            ! Error flag (0=success, -1=invalid input)
  
  integer :: row_sizes(N)                 ! Number of elements in each row (1-based)
  integer :: prefix_sums(N)               ! Cumulative sum of elements
  real :: target, target_cumulative       ! Target elements per process
  integer :: total_elements               ! Total elements in the upper triangle
  integer :: current_start, proc_id, j, best_split, min_diff, current_sum, p

  ! Initialize error flag
  ierr = 0

  ! Validate input
  if (N < 1 .or. num_processes < 1) then
    ierr = -1
    return
  end if

  ! Compute row sizes (1-based)
  do j = 1, N
    row_sizes(j) = N - j + 1
  end do

  ! Compute prefix sums (cumulative elements up to each row)
  prefix_sums(1) = row_sizes(1)
  do j = 2, N
    prefix_sums(j) = prefix_sums(j-1) + row_sizes(j)
  end do
  total_elements = prefix_sums(N)

  ! Calculate target elements per process
  target = real(total_elements) / real(num_processes)

  current_start = 1  ! Start from the first row

  ! Distribute rows across processes
  do proc_id = 1, num_processes
    if (proc_id == num_processes) then
      ! Last process gets all remaining rows
      start_rows(proc_id) = current_start
      end_rows(proc_id) = N
      exit
    end if

    ! Target cumulative elements for this process
    target_cumulative = proc_id * target

    ! Find the best split point
    min_diff = HUGE(0)  ! Initialize to a large number
    best_split = current_start

    do j = current_start, N
      current_sum = prefix_sums(j)
      if (ABS(current_sum - target_cumulative) < min_diff) then
        min_diff = ABS(current_sum - target_cumulative)
        best_split = j
      else if (current_sum > target_cumulative) then
        exit  ! No better split exists beyond this point
      end if
    end do

    ! Assign rows to the current process
    start_rows(proc_id) = current_start
    end_rows(proc_id) = best_split
    current_start = best_split + 1

    ! Exit early if no more rows to assign (handle remaining processes)
    if (current_start > N) then
      ! Set remaining processes to have no rows
      do p = proc_id + 1, num_processes
        start_rows(p) = 1
        end_rows(p) = 0
      end do
      exit  ! Break outer loop after handling remaining processes
    end if
  end do

  ! Calculate elements per process
  do proc_id = 1, num_processes
    if (start_rows(proc_id) > end_rows(proc_id)) then
      elements(proc_id) = 0
    else
      elements(proc_id) = sum(row_sizes(start_rows(proc_id):end_rows(proc_id)))
    end if
  end do

  ! Verify all elements are assigned
  if (sum(elements) /= total_elements) ierr = -2

end subroutine distribute_rows


program test_distribute_rows
    implicit none
    integer, parameter :: N = 20, num_processes = 3
    integer :: start_rows(num_processes), end_rows(num_processes), elements(num_processes)
    integer :: ierr, i

    call distribute_rows(N, num_processes, start_rows, end_rows, elements, ierr)

    if (ierr /= 0) then
      print *, "Error code:", ierr
      stop
    end if

    do i = 1, num_processes
      if (start_rows(i) > end_rows(i)) then
        print "('Process ', I0, ': No rows assigned')", i
      else
        print "('Process ', I0, ': rows ', I0, '-', I0, ', elements=', I0)", &
              i, start_rows(i), end_rows(i), elements(i)
      end if
    end do
end program test_distribute_rows