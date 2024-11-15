include "mkl_spblas.f90"
! Sparse BLAS matrix-vector product example
!
! The following code demonstrates the creation of
! a CSR matrix, and the sparse matrix-vector product
! using the sparse BLAS functions in Intel MKL.
!
! For linking options for your specific platform see
! the Intel MKL link line advisor located at:
!   https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html
!
! In my case I was able to compile the example with the 
! following code, by executing the following commands
! from a shell script:
!
! LINK="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl"
! INC=-I"${MKLROOT}/include/"
! ifort -warn all -O3 $INC/mkl_spblas.f90 test_spblas.f90 -o test_spblas $LINK
!
! If you have the MKL pkg-config configured then a command like
! the one below should do:
!
! $ ifort test_spblas.f90 -o test_spblas `pkg-config --libs mkl-static-lp64`
!
program test_spblas

    use mkl_spblas
    implicit none
  
  
    integer, parameter :: rows = 4
    integer, parameter :: cols = 6
  
    integer, parameter :: nnz = 8
  
    integer :: ia(rows+1), ja(nnz), stat
    real :: values(nnz), x(6), y(4)
  
    type(sparse_matrix_t) :: a
    type(matrix_descr) :: descr
  
  
    ! Matrix example taken from: 
    ! https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
    !
    !     | 10  20  0  0  0  0 |
    ! A = |  0  30  0 40  0  0 |
    !     |  0   0 50 60 70  0 |
    !     |  0   0  0  0  0 80 | 
   
    ia = [1,3,5,8,9]
    ja = [1,2,2,4,3,4,5,6]
    values = [10, 20, 30, 40, 50, 60, 70, 80]
  
    stat = mkl_sparse_s_create_csr(a,SPARSE_INDEX_BASE_ONE,rows,cols,ia(1:4),ia(2:5),ja,values)
    print *, "stat create = ", stat
  
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
  
    x = [1,1,1,1,1,1]
    stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.0,a,descr,x,0.0,y)
    print *, "stat mv = ", stat
  
    print *, "result   = ", y
    print *, "expected = ", [30., 70., 180., 80.]
  
  end program