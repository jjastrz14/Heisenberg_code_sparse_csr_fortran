!test for mm_csr in BLAS
subroutine mmm_csr_test()

    use mkl_spblas
    use iso_C_binding
    implicit none 

    integer :: i, j, stat_permutation, stat_H_csr, stat_csr_mutliplication, info, stat_export, nrow, ncol, info_spmmd
    double precision, allocatable :: values_array_H(:), values_array_P(:)
    integer, allocatable :: ia_H_B(:), ia_H_E(:), ja_H(:), ia_P_B(:), ia_P_E(:), ja_P(:)

    !BLAS types
    type(sparse_matrix_t) :: p_matrix_csr
    type(sparse_matrix_t) :: H_matrix_csr
    type(sparse_matrix_t) :: H_matrix_permuted
    type(sparse_matrix_t) :: H_permuted
    type(matrix_descr) :: descr_H, descr_P
    !C_pointer types
    integer(C_INT) :: indexing
    type(C_PTR)   :: ia_start, ia_end, ja, values

    integer, pointer :: ia_export_B(:), ia_export_E(:), ja_export(:)
    double precision, pointer :: values_array_export(:)

    write(*,*) '----------- Sparse BLAS 4x4 matrix test:-----------'
    write(*,*) " "
    allocate(values_array_P(4), ia_P_B(4), ia_P_E(4), ja_P(4))

    ! P matrix
    ! 1 0 0 0
    ! 0 0 1 0
    ! 0 1 0 0
    ! 0 0 0 1 
    !values vector
    values_array_P(1) = 1.0d0
    values_array_P(2) = 1.0d0
    values_array_P(3) = 1.0d0
    values_array_P(4) = 1.0d0
    !columns vector
    ja_P(1) = 1
    ja_P(2) = 3
    ja_P(3) = 2
    ja_P(4) = 4
    ! csr pointer B
    ia_P_B(1) = 1;
    ia_P_B(2) = 2;
    ia_P_B(3) = 3;
    ia_P_B(4) = 4;
    ! csr pointer E
    ia_P_E(1) = 2;
    ia_P_E(2) = 3;
    ia_P_E(3) = 4;
    ia_P_E(4) = 5;


    ! H matrix symmetric upper triangler to csr3
    allocate(values_array_H(6), ia_H_B(4), ia_H_E(4), ja_H(6))
    ! H matrix
    ! 1 0 1 0
    ! 0 1 0 2
    ! 1 0 1 0
    ! 0 2 0 1
    !values vector
    values_array_H(1) = 1.0d0
    values_array_H(2) = 1.0d0
    values_array_H(3) = 1.0d0
    values_array_H(4) = 2.0d0
    values_array_H(5) = 1.0d0
    values_array_H(6) = 1.0d0
    !columns vector
    ja_H(1) = 1
    ja_H(2) = 3
    ja_H(3) = 2
    ja_H(4) = 4
    ja_H(5) = 3
    ja_H(6) = 4
    ! csr pointer B
    ia_H_B(1) = 1;
    ia_H_B(2) = 3;
    ia_H_B(3) = 5;
    ia_H_B(4) = 6;
    ! csr pointer E
    ia_H_E(1) = 3;
    ia_H_E(2) = 5;
    ia_H_E(3) = 6;
    ia_H_E(4) = 7;

    !rewriting permutation matrix to intel mkl internal csr format
    stat_permutation = mkl_sparse_d_create_csr(p_matrix_csr, SPARSE_INDEX_BASE_ONE, 4, 4 , ia_P_B, ia_P_E, ja_P, values_array_P)
    print *, "stat CSR Permutation matrix create = ", stat_permutation

    stat_export = mkl_sparse_d_export_csr(p_matrix_csr, indexing, nrow, ncol, ia_start, ia_end, ja, values)
    print *, "Export P = ", stat_export

    !   Converting C into Fortran pointers
    call C_F_POINTER(ia_start, ia_export_B, [nrow])
    call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
    call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
    call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
    
    write(*,*)  ' '
    write(*,*) 'Exported CSR format of P: '
    print *, 'values array: ', values_array_export
    print *, 'ja_P: ', ja_export
    print *, 'ia_P; pointerB: ', ia_export_B
    print *, 'ia_P; pointerE: ', ia_export_E

    descr_P%type = SPARSE_MATRIX_TYPE_GENERAL
    info = MKL_SPARSE_OPTIMIZE(p_matrix_csr)
    write(*,*) 'info optimize P', info

    !rewriting H matrix to intel mkl internal csr format
    stat_H_csr = mkl_sparse_d_create_csr(H_matrix_csr, SPARSE_INDEX_BASE_ONE, 4, 4 , ia_H_B, ia_H_E, ja_H, values_array_H)
    print *, "stat CSR H matrix create = ", stat_H_csr

    descr_H % type = SPARSE_MATRIX_TYPE_SYMMETRIC ! it describes H matrix, P is general (full csr) by assumption of mkl_sparse_sypr
    descr_H % mode = SPARSE_FILL_MODE_UPPER
    descr_H % diag = SPARSE_DIAG_NON_UNIT

    info = MKL_SPARSE_OPTIMIZE(H_matrix_csr)
    write(*,*) 'info optimize H', info
    !info = mkl_sparse_order(H_matrix_csr)
    !write(*,*) 'info optimize P', info

    stat_export = mkl_sparse_d_export_csr(H_matrix_csr, indexing, nrow, ncol, ia_start, ia_end, ja, values)
    print *, "Export H_csr = ", stat_export

    !   Converting C into Fortran pointers
    call C_F_POINTER(ia_start, ia_export_B, [nrow])
    call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
    call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
    call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
    
    write(*,*)  ' '
    write(*,*) 'Exported CSR format of H not permuted: '
    print *, 'values array: ', values_array_export
    print *, 'ja_H: ', ja_export
    print *, 'ia_H; pointerB: ', ia_export_B
    print *, 'ia_H; pointerE: ', ia_export_E


    !Doing a product of three csr matrices P @ H @ P.T <-SPARSE_OPERATION_NON_TRANSPOSE
    stat_csr_mutliplication = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, p_matrix_csr, H_matrix_csr, descr_H, H_matrix_permuted, SPARSE_STAGE_FULL_MULT)

    print *, "stat P @ H matrix @ P.T  = ", stat_csr_mutliplication

    info = mkl_sparse_order(H_matrix_permuted)
    stat_export = mkl_sparse_d_export_csr(H_matrix_permuted, indexing, nrow, ncol, ia_start, ia_end, ja, values)

    print *, "Export H permuted  = ", stat_export
    
    ! H permuted
    ! 1 1 0 0
    ! 1 1 0 0
    ! 0 0 1 2
    ! 0 0 2 1

    !   Converting C into Fortran pointers
    call C_F_POINTER(ia_start, ia_export_B, [nrow])
    call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
    call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
    call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
    
    write(*,*)  ' '
    write(*,*) 'Exported CSR format of H permuted: '
    print *, 'values array: ', values_array_export
    print *, 'ja_H_permuted: ', ja_export
    print *, 'ia_H_permuted; pointerB: ', ia_export_B
    print *, 'ia_H_permuted; pointerE: ', ia_export_E
    
    !   Release internal representation of CSR matrix
    info = mkl_sparse_destroy(p_matrix_csr)
    info = mkl_sparse_destroy(H_matrix_csr)
    info = mkl_sparse_destroy(H_matrix_permuted)

    deallocate(values_array_P, ia_P_B, ia_P_E, ja_P)
    deallocate(values_array_H, ia_H_B, ia_H_E, ja_H)

    write(*,*) '----------- END Sparse BLAS test: -----------'
    write(*,*) " "

end subroutine mmm_csr_test