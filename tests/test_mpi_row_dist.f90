subroutine distribute_rows(N, num_processes, start_rows, end_rows, elements, ierr)
  implicit none
  integer(8), intent(in) :: N                ! Size of the matrix (N x N)
  integer(8), intent(in) :: num_processes    ! Number of MPI processes
  integer(8), intent(out) :: start_rows(num_processes) ! Start row (1-based)
  integer(8), intent(out) :: end_rows(num_processes)   ! End row (1-based)
  integer(8), intent(out) :: elements(num_processes)   ! Elements per process
  integer, intent(out) :: ierr            ! Error flag (0=success, -1=invalid input)

  double precision :: total_elements               ! Total elements in upper triangle
  double precision :: target                           ! Target elements per process
  integer(8) :: current_start, proc_id, j, best_split, min_diff, temp_sum, p
  integer(8) :: n_rows, first_element, last_element

  ! Initialize error flag
  ierr = 0

  ! Validate input
  if (N < 1 .or. num_processes < 1) then
    ierr = -1
    return
  end if

  ! Total elements = N*(N+1)/2
  total_elements = N * (N + 1.0d0) / 2.0d0
  target = (total_elements+0.0d0) / (num_processes+0.0d0)
  write(*,*) "Total elements:", total_elements
  write(*,*) "Target elements per process:", target

  current_start = 1  ! Start from row 1

  ! Distribute rows across processes
  do proc_id = 1, num_processes
    if (proc_id == num_processes) then
      ! Last process gets remaining rows
      start_rows(proc_id) = current_start
      end_rows(proc_id) = N
      exit
    end if

    min_diff = HUGE(0)  ! Initialize to largest integer
    best_split = current_start
    temp_sum = 0

    ! Dynamically compute cumulative sum
    do j = current_start, N
      temp_sum = temp_sum + (N - j + 1)  ! Elements in row j
      if (abs(temp_sum - target) < min_diff) then
        min_diff = abs(temp_sum - target)
        best_split = j
      else if (temp_sum > target) then
        exit  ! Stop if sum exceeds target
      end if
    end do

    ! Assign rows to the current process
    start_rows(proc_id) = current_start
    end_rows(proc_id) = best_split
    current_start = best_split + 1

    ! Handle no remaining rows
    if (current_start > N) then
      do p = proc_id + 1, num_processes
        start_rows(p) = 1
        end_rows(p) = 0
      end do
      exit
    end if
  end do

  ! Calculate elements using arithmetic series formula
  do proc_id = 1, num_processes
    if (start_rows(proc_id) > end_rows(proc_id)) then
      elements(proc_id) = 0
    else
      n_rows = end_rows(proc_id) - start_rows(proc_id) + 1
      first_element = N - start_rows(proc_id) + 1
      last_element = N - end_rows(proc_id) + 1
      elements(proc_id) = n_rows * (first_element + last_element) / 2
    end if
  end do

  ! Verify total elements
  if (sum(elements) /= total_elements) ierr = -2

end subroutine distribute_rows


program test_distribute_rows
    implicit none
    integer(8), parameter :: N = 184756 , num_processes = 3 !252 !184756
    integer(8) :: start_rows(num_processes), end_rows(num_processes), elements(num_processes)
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