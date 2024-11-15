include "mkl_vsl.f90" !necessary for intel mkl random numbers
 
module omp_mkl_test
    contains
    subroutine omp_mkl_small_test()
        use omp_lib
        use mkl_vsl
        implicit none

        integer :: i, j, omp_id, brng, rand_method, rand_seed, errcode, no_of_rand_num
        real(kind=8) rand_start, rand_end
        real(kind=8), allocatable :: rand_num(:)
        type (VSL_STREAM_STATE) :: stream(1)

        CHARACTER*1 :: jobz, range, uplo
        integer :: N_mat, il, iu, ldz, liwork, m, lwork, lrwork, info
        integer, allocatable :: iwork(:), isuppz(:)
        double precision :: vl, vu, abstol
        double complex, allocatable :: H_test(:,:), work_complex(:), eigen_vec(:,:)
        double precision, allocatable :: rwork(:), w(:)
        double complex :: norm

        double complex, allocatable  :: inv_H_test(:,:)
        integer, allocatable :: ipv(:)

        !   0 - part checking OpenMP + IntelMKL (oneAPI) random numbers/diag/inversion/
        write(*,*) 'OpenMP processors check'
        !$OMP PARALLEL PRIVATE(omp_id)
        omp_id = omp_get_thread_num()
        write(*,*) 'OpenMP thread check', omp_id
        !$OMP END PARALLEL

    

        write(*,*) 'Math Kernel Library random numbers check'
        ! *** random number generator parameters *** !
        brng   = VSL_BRNG_MT2203
        rand_method = VSL_RNG_METHOD_UNIFORM_STD
        rand_start = 0.0
        rand_end = 1.0
        rand_seed = 1
        no_of_rand_num = 10
        allocate(rand_num(no_of_rand_num))
        errcode = vslnewstream(stream(1), brng, rand_seed)
        write(*,*) 'random stream init code: ', errcode
        errcode=vdrnguniform( rand_method, stream(1), no_of_rand_num, rand_num, rand_start, rand_end )
        write(*,*) 'random numbers error code: ', errcode
        do i = 1, no_of_rand_num
            write(*,*) i, rand_num(i)
        end do
        deallocate(rand_num)
        errcode=vsldeletestream(stream(1))
        write(*,*) 'random stream kill code: ', errcode
        !end do
        write(*,*) 'Small Hermitean matrix diagonalization test'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        N_mat = 3 ! The order of the matrix
        vl = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il=1  !the index of the smallest eigenvalue to be returned.
        iu=N_mat !the index of the largest eigenvalue to be returned.
        abstol=10**(-10.0d0)
        ldz=N_mat
        liwork= 10*N_mat
        lwork=2*N_mat
        lrwork=24*N_mat
        allocate (H_test(N_mat,N_mat))
        H_test(1,1) = dcmplx(1.0d0, 0.0d0)
        H_test(1,2) = dcmplx(1.0d0, 1.0d0)
        H_test(1,3) = dcmplx(2.0d0, 2.0d0)
        H_test(2,1) = dcmplx(1.0d0,-1.0d0)
        H_test(2,2) = dcmplx(2.0d0, 0.0d0)
        H_test(2,3) = dcmplx(3.0d0, 3.0d0)
        H_test(3,1) = dcmplx(2.0d0,-2.0d0)
        H_test(3,2) = dcmplx(3.0d0,-3.0d0)
        H_test(3,3) = dcmplx(3.0d0, 0.0d0)
        allocate ( eigen_vec(N_mat,N_mat), work_complex(lwork),iwork(liwork),w(N_mat),isuppz(2*N_mat), rwork(lrwork)  )
        call zheevr(jobz, range, uplo, N_mat, H_test, N_mat, vl, vu, il, iu, abstol, m, w, eigen_vec, ldz, isuppz, work_complex, lwork, rwork, lrwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work_complex(1)
        write(*,*) "number of eigeval found:", m
        do i= 1, m
            write(*,*) i, 'th eigenval is', w(i)
        end do
        write(*,*) 'compare with Matlab: -2.4801476299237900'
        write(*,*) 'compare with Matlab:  0.5056116166572769'
        write(*,*) 'compare with Matlab:  7.9745360132665200'

    

        write(*,*) 'eigenvec norm test:'
        do i=1,m
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1,N_mat
                norm = norm + eigen_vec(j,i)*conjg(eigen_vec(j,i))
            end do
            write(*,*) i, norm
        end do
        deallocate ( eigen_vec, work_complex, iwork, w, isuppz, rwork  )
        !
        write(*,*) 'Small Hermitean matrix inversion test'
        H_test(1,1) = dcmplx(1.0d0, 0.0d0)
        H_test(1,2) = dcmplx(1.0d0, 1.0d0)
        H_test(1,3) = dcmplx(2.0d0, 2.0d0)
        H_test(2,1) = dcmplx(1.0d0,-1.0d0)
        H_test(2,2) = dcmplx(2.0d0, 0.0d0)
        H_test(2,3) = dcmplx(3.0d0, 3.0d0)
        H_test(3,1) = dcmplx(2.0d0,-2.0d0)
        H_test(3,2) = dcmplx(3.0d0,-3.0d0)
        H_test(3,3) = dcmplx(3.0d0, 0.0d0)
        lwork = N_mat
        allocate(inv_H_test(N_mat,N_mat), ipv(N_mat), work_complex(lwork) )
        call zgetrf(N_mat, N_mat, H_test, N_mat, ipv, info ) !number of operations: (8/3)*N_mat^3
        write(*,*) 'blas LU info ', info
        call zgetri(N_mat, H_test, N_mat, ipv, work_complex, N_mat, info )
        !number of operations: (16/3)*M^3
        write(*,*) 'blas inversion info ', info
        write(*,*) lwork, 'inversion optimal lwork ', work_complex(1)
        do i = 1, N_mat
            do j = 1, N_mat
                write(*,*) i, j, H_test(i,j)
            end do
        end do
        write(*,*) 'inverted matrix should be:'
        write(*,*) ' 1.2+0.0i, -0.9+0.3i,  0.4-0.2i'
        write(*,*) '-0.9-0.3i,  0.5+0.0i, -0.1+0.3i'
        write(*,*) ' 0.4+0.2i, -0.1-0.3i,  0.0+0.0i'
        deallocate(H_test, ipv, work_complex)

    end subroutine omp_mkl_small_test

end module omp_mkl_test

module spin_systems
    contains

    subroutine str(k, name)
        implicit none 

        !   "Convert an integer to string."

        character(len=20), intent(out):: name
        integer, intent(in) :: k 

        write(name, *) k 
        name = adjustl(name)

    end subroutine str

    subroutine H_XXX_full()
        use omp_lib
        use mkl_vsl
        use omp_mkl_test
        implicit none
        
        integer :: i, j, N_spin, J_spin, N_spin_max
        double complex, allocatable :: H_full(:,:), Sp_site_ind_Sm_site_ind_p1(:,:), Sm_site_ind_Sp_site_ind_p1(:,:), Sz_site_ind_Sz_site_ind_p1(:,:)  
        integer :: site_ind, m, n, p, q, r, s
        integer :: ind_i, ind_j, a_i, a_j, b_i, b_j, c_i, c_j
        double precision :: A_val, B_val, C_val, result_ij, norm
 
        CHARACTER*1 :: jobz, range, uplo
        integer :: il, iu, ldz, liwork, lwork, lrwork, info, m_eig
        double precision :: vl, vu, abstol
        double complex, allocatable :: work_complex(:), eigen_vec(:,:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: rwork(:), w(:)
        
        open (unit=10, file="H_check.dat", recl=512)
        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        
        call omp_mkl_small_test()

        N_spin = 4
        J_spin = 1
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
    
    
        allocate( H_full(N_spin_max, N_spin_max), Sp_site_ind_Sm_site_ind_p1(N_spin_max, N_spin_max), Sm_site_ind_Sp_site_ind_p1(N_spin_max, N_spin_max), Sz_site_ind_Sz_site_ind_p1(N_spin_max, N_spin_max) )
    
        H_full = dcmplx(0.0d0,0.0d0)
       
        write(*,*) 'H_full(1,1)', H_full(1,1)

    
        do site_ind = 1, (N_spin-1)
            m = 2**(site_ind-1)
            n = m
            p = 4 ! N_spin independent because this referes to S+S- operator for spin 1/2 
            q = 4
            r = 2**(N_spin-(site_ind+1))
            s = r
            !write(*,*) site_ind, m, n, p, q, r, s
        
            Sp_site_ind_Sm_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            Sm_site_ind_Sp_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            Sz_site_ind_Sz_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            
            do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
                !write(*,*) 'site_ind, ind_i, ind_j = ',  site_ind, ind_i, ind_j
                !i = ind_i-1
                !j = ind_j-1
                a_i = FLOOR( (ind_i-1.0d0)/( (p*r)+0.0d0 ) )
                a_j = FLOOR( (ind_j-1.0d0)/( (q*s)+0.0d0 ) )
                b_i = MOD(FLOOR( (ind_i-1.0d0)/( r+0.0d0 ) ), p)
                b_j = MOD(FLOOR( (ind_j-1.0d0)/( s+0.0d0 ) ), q)
                c_i = MOD(ind_i-1, r)
                c_j = MOD(ind_j-1, s)
                                
                A_val = 0.0d0
                B_val = 0.0d0
                C_val = 0.0d0
            
                if (a_i == a_j) then
                    A_val = 1.0d0
                else
                    A_val = 0.0d0
                end if
            
                if (c_i == c_j) then
                    C_val = 1.0d0
                else
                    C_val = 0.0d0
                end if
            
                ! For S+S-:
                if (b_i == 1 .AND. b_j == 2) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if
                result_ij = A_val*B_val*C_val
                Sp_site_ind_Sm_site_ind_p1(ind_i,ind_j) = result_ij
            
                ! For S-S+:
                if (b_i == 2 .AND. b_j == 1) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if
                result_ij = A_val*B_val*C_val
                Sm_site_ind_Sp_site_ind_p1(ind_i,ind_j) = result_ij
            
                ! For SzSz:
                if (b_i == 0 .AND. b_j == 0) then
                    B_val = 0.250d0
                else if (b_i == 1 .AND. b_j == 1) then
                    B_val = -0.250d0
                else if (b_i == 2 .AND. b_j == 2) then
                    B_val = -0.250d0
                else if (b_i == 3 .AND. b_j == 3) then
                    B_val = 0.250d0
                else
                    B_val = 0.0d0
                end if
                result_ij = A_val*B_val*C_val
                Sz_site_ind_Sz_site_ind_p1(ind_i,ind_j) = result_ij
            
            end do
            end do
        
            H_full = H_full + 1.0d0/2.0d0 *(Sp_site_ind_Sm_site_ind_p1+Sm_site_ind_Sp_site_ind_p1) + Sz_site_ind_Sz_site_ind_p1
            
        end do
        
    
        write(*,*) 'H_full matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il=1  !the index of the smallest eigenvalue to be returned.
        iu=N_spin_max !the index of the largest eigenvalue to be returned.
        abstol=10**(-10.0d0)
        ldz=N_spin_max
        liwork= 10*N_spin_max
        lwork = 2*N_spin_max !check with work complex(1)
        lrwork= 24*N_spin_max

        allocate (eigen_vec(N_spin_max,N_spin_max), work_complex(lwork),iwork(liwork),w(N_spin_max),isuppz(2*N_spin_max), rwork(lrwork))
        write(*,*) 'just before allocation'
        call zheevr(jobz, range, uplo, N_spin_max, H_full, N_spin_max, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, ldz, isuppz, work_complex, lwork, rwork, lrwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work_complex(1)
        write(*,*) "number of eigeval found:", m_eig
        do i= 1, m_eig
            write(*,*) i, 'th eigenval is', w(i)
        end do
    
        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, N_spin_max
                norm = norm + eigen_vec(j,i)*conjg(eigen_vec(j,i))
            end do
            write(*,*) i, norm
        end do
        deallocate ( eigen_vec, work_complex, iwork, w, isuppz, rwork  )
    
        deallocate( H_full, Sp_site_ind_Sm_site_ind_p1, Sm_site_ind_Sp_site_ind_p1, Sz_site_ind_Sz_site_ind_p1 )
        
        close(10)


    end subroutine H_XXX_full

    subroutine H_XXX_filling(N_spin, J_spin, i, j, H_full_ij)
        implicit none
        integer, intent(in) :: N_spin, i, j
        double precision, intent(in) :: J_spin
        double precision, intent(out) :: H_full_ij
        integer :: N_spin_max  
        double precision :: Sp_site_ind_Sm_site_ind_p1_ij, Sm_site_ind_Sp_site_ind_p1_ij, Sz_site_ind_Sz_site_ind_p1_ij
        integer :: site_ind, m, n, p, q, r, s
        integer :: ind_i, ind_j, a_i, a_j, b_i, b_j, c_i, c_j
        double precision :: A_val, B_val, C_val

        ! N spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        N_spin_max = 2**N_spin

        H_full_ij = 0.0d0
                           
        do site_ind = 1, (N_spin-1)
            m = 2**(site_ind-1)
            n = m
            p = 4 ! N_spin independent because this referes to S+S- operator for spin 1/2 
            q = 4
            r = 2**(N_spin-(site_ind+1))
            s = r
            !write(*,*) site_ind, m, n, p, q, r, s
        
            Sp_site_ind_Sm_site_ind_p1_ij = 0.0d0
            Sm_site_ind_Sp_site_ind_p1_ij = 0.0d0
            Sz_site_ind_Sz_site_ind_p1_ij = 0.0d0
            
            !do ind_i = 1, N_spin_max
            !do ind_j = 1, N_spin_max
            ind_i = i
            ind_j = j
                !write(*,*) 'site_ind, ind_i, ind_j = ',  site_ind, ind_i, ind_j
                !i = ind_i-1
                !j = ind_j-1
                a_i = FLOOR( (ind_i-1.0d0)/( (p*r)+0.0d0 ) )
                a_j = FLOOR( (ind_j-1.0d0)/( (q*s)+0.0d0 ) )
                b_i = MOD(FLOOR( (ind_i-1.0d0)/( r+0.0d0 ) ), p)
                b_j = MOD(FLOOR( (ind_j-1.0d0)/( s+0.0d0 ) ), q)
                c_i = MOD(ind_i-1, r)
                c_j = MOD(ind_j-1, s)
                                
                A_val = 0.0d0
                B_val = 0.0d0
                C_val = 0.0d0
            
                if (a_i == a_j) then
                    A_val = 1.0d0
                else
                    A_val = 0.0d0
                end if
            
                if (c_i == c_j) then
                    C_val = 1.0d0
                else
                    C_val = 0.0d0
                end if
            
                ! For S+S-:
                if (b_i == 1 .AND. b_j == 2) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if
                Sp_site_ind_Sm_site_ind_p1_ij = A_val*B_val*C_val
            
                ! For S-S+:
                if (b_i == 2 .AND. b_j == 1) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if                    
                Sm_site_ind_Sp_site_ind_p1_ij = A_val*B_val*C_val
            
                ! For SzSz:
                if (b_i == 0 .AND. b_j == 0) then
                    B_val = 0.250d0
                else if (b_i == 1 .AND. b_j == 1) then
                    B_val = -0.250d0
                else if (b_i == 2 .AND. b_j == 2) then
                    B_val = -0.250d0
                else if (b_i == 3 .AND. b_j == 3) then
                    B_val = 0.250d0
                else
                    B_val = 0.0d0
                end if

                Sz_site_ind_Sz_site_ind_p1_ij= A_val*B_val*C_val
            
            !end do
            !end do
        
            H_full_ij = H_full_ij + 1.0d0/2.0d0 *(Sp_site_ind_Sm_site_ind_p1_ij+Sm_site_ind_Sp_site_ind_p1_ij) + Sz_site_ind_Sz_site_ind_p1_ij
            
        end do


    end subroutine H_XXX_filling

    
    subroutine H_XXX_diag(N_spin, J_spin)
        use omp_lib
        use mkl_vsl
        use omp_mkl_test
        implicit none

        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer :: i, j, ind_i, ind_j, N_spin_max
        double precision, allocatable :: H_full(:,:)
        double precision :: norm, H_full_ij, J_spin_number
 
        CHARACTER*1 :: jobz, range, uplo
        integer :: il, iu, ldz, liwork, lwork, info, m_eig
        double precision :: vl, vu, abstol
        double precision, allocatable :: work(:), eigen_vec(:,:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: w(:)
        
        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        ! Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

        file_name1 = 'H_eigenvals_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name2 = 'H_eigenvectors_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name3 = 'H_norm_test_ful_' // trim(adjustl(N_spin_charachter)) // '.dat'

        open (unit=10, file= trim(file_name1), recl=512)
        open (unit=11, file=trim(file_name2), recl=512)
        open (unit=12, file= trim(file_name3), recl=512)
        
        ! open (unit=11, file="H_full.dat", recl=512)
        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
    
        call omp_mkl_small_test()
        
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        !call H_XXX_filling(N_spin, J_spin, i, j, H_full_ij)

        allocate( H_full(N_spin_max, N_spin_max))
        J_spin_number = J_spin
        H_full = 0.0d0
       
        !!! H_full filling
        write(*,*) 'H_full: '
        do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
                call H_XXX_filling(N_spin, J_spin_number, ind_i, ind_j, H_full_ij)
                H_full(ind_i,ind_j) = H_full_ij
                ! write(11,*) ind_i, ind_j, H_full(ind_i,ind_j)
            end do
        end do
        
    
        write(*,*) 'H_full matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl     = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu     = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il     = 1  !the index of the smallest eigenvalue to be returned.
        iu     = N_spin_max !the index of the largest eigenvalue to be returned.
        abstol = 10**(-10.0d0)
        ldz    = N_spin_max
        lwork  = 26*N_spin_max
        liwork = 10*N_spin_max

        allocate ( w(N_spin_max), eigen_vec(N_spin_max,N_spin_max), isuppz(2*N_spin_max), work(lwork), iwork(liwork) )
        write(*,*) 'just before allocation'
       
        call dsyevr(jobz, range, uplo, N_spin_max, H_full, N_spin_max, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, N_spin_max, isuppz, work, lwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work(1)
        write(*,*) "number of eigeval found:", m_eig

        write(10,*) 'i-th , eigenvalue'
        do i= 1, m_eig
            write(10,*) i, ',', w(i) 
        end do

        write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        do i=1, m_eig
            do j=1, N_spin_max
                write(11,*) eigen_vec(j,i)
            end do
            write(11,*) " "
        end do

        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, N_spin_max
                norm = norm + eigen_vec(j,i)*(eigen_vec(j,i))
            end do
            write(12,*) i, norm
        end do
        
        deallocate (H_full,  w, eigen_vec, isuppz, work, iwork )
                            
        close(10)
        close(11)
        close(12)

    end subroutine H_XXX_diag

    subroutine sparse_dfeast_test()
        
        implicit none
        double complex, allocatable :: test_matrix(:,:)
        CHARACTER*1 :: jobz, uplo            
        integer :: i, n, lda, lwork, info
        double precision, allocatable :: eigenvalues(:), rwork(:)
        double complex, allocatable ::  work(:,:)
        double precision, allocatable :: values_array(:), x(:,:)
        integer, allocatable :: ia(:), ja(:)
        integer :: fpm(128)
        double precision :: emin, emax
        integer :: m0
        double precision :: epsout
        integer :: loop
        double precision, allocatable :: e(:)
        integer :: m
        double precision, allocatable :: res(:)
        
        write(*,*) 'FEAST testing'
        n = 5
        allocate (test_matrix(n,n))
        test_matrix = dcmplx(0.0d0, 0.0d0)
        test_matrix(1,1) = dcmplx(1.0d0, 0.0d0)
        test_matrix(2,2) = dcmplx(0.0d0, 0.0d0)
        test_matrix(3,3) = dcmplx(2.0d0, 0.0d0)
        test_matrix(4,4) = dcmplx(0.0d0, 0.0d0)
        test_matrix(5,5) = dcmplx(3.0d0, 0.0d0)
        test_matrix(1,4) = dcmplx(5.0d0, 0.0d0)
        test_matrix(2,3) = dcmplx(6.0d0, 0.0d0)
        test_matrix(3,5) = dcmplx(7.0d0, 0.0d0)
        
        !standard diagonalziation
        jobz='V' !N - only eigenvalues, V - eigenvalues + eigenvectors
        uplo='U' !upper hermitean matrix part passed to routine
        n=5
        lda=n
        lwork=2*n-1  
        allocate(work(lwork,lwork),rwork(3*n-2),eigenvalues(n))
        
        call zheev(jobz, uplo, n, test_matrix, lda, eigenvalues, work, lwork, rwork, info)
        if (info /= 0) then
            write(*,*) 'problem with zheev, info=', info
        end if 
            
        do i=1,n
            write(*,*) i, eigenvalues(i) 
        enddo     
        deallocate(work,rwork,eigenvalues)

        write(*,*) 'Zheev diagonalization finished, FEAST start'

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        allocate( values_array(8), ia(6), ja(8) )
         
        values_array(1) = dcmplx(1.0d0, 0.0d0)
        values_array(2) = dcmplx(5.0d0, 0.0d0)
        values_array(3) = dcmplx(0.0d0, 0.0d0)
        values_array(4) = dcmplx(6.0d0, 0.0d0)
        values_array(5) = dcmplx(2.0d0, 0.0d0)
        values_array(6) = dcmplx(7.0d0, 0.0d0)
        values_array(7) = dcmplx(0.0d0, 0.0d0)
        values_array(8) = dcmplx(3.0d0, 0.0d0)
        
        ja(1) = 1
        ja(2) = 4
        ja(3) = 2
        ja(4) = 3
        ja(5) = 3
        ja(6) = 5
        ja(7) = 4
        ja(8) = 5
        
        ia(1) = 1
        ia(2) = 3
        ia(3) = 5
        ia(4) = 7
        ia(5) = 8
        ia(6) = 9 !no of non-zero +1 !!!
           
        
        !n=n     ! Sets the size of the problem
        !a=non_zero_array     ! Array containing the nonzero elements of the upper triangular part of the matrix A
        call feastinit (fpm)  ! function specifying default parameters fpm of FEAST algorithm
        fpm(1) = 1
        fpm(2) = 12 !can be more, can be less
        fpm(3) = 8 !eps 10^-fpm(3)
        fpm(4) = 20 ! max number of feast loops
        fpm(5) = 0 !initial subspace
        fpm(6) = 0! stopping criterion
        fpm(7) = 5 !Error trace sigle prec stop crit
        fpm(14) = 0 ! standard use of feast
        fpm(27) = 1 !check input matrices
        fpm(28) = 1 !check if B is positive definite?         
        uplo='U' ! If uplo = 'U', a stores the upper triangular parts of A.
        emin = -20.0d0 ! The lower ... &
        emax =  20.0d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = 5 !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
           !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
           !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = 5
        allocate( x(n,m0), e(m0), res(m0) )
        write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_hcsrev... '
        call dfeast_scsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) 'dfeast_hcsrev info=', info
        write(*,*) 'After dfeast_hcsrev... '
            
        if (info /= 0) then
           write(*,*) 'problem with dfeast_hcsrev, info=', info
        end if 
        
        write(*,*) 'no of eigenvalues found= ', m
        do i = 1, m
            write(*,*) i, e(i)
        end do
                    
        deallocate(test_matrix, values_array, ia, ja, x, e, res)

    end subroutine sparse_dfeast_test


    subroutine H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero)
        use omp_lib
        use mkl_vsl
        use omp_mkl_test
        implicit none
        
        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer, intent(out) :: no_of_nonzero
        integer :: i, j, ind_i, ind_j,  N_spin_max, counter, ind_ia            
        double precision :: norm, H_full_ij

        open (unit=20, file="ia.dat", recl=512)
        open (unit=21, file="ja.dat", recl=512)
        open (unit=22, file="values_array.dat", recl=512)
  
        ! N spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1

        !N_spin = 4
        !J_spin = 1
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        ! we start with saving ia, ja and val_arr into files (needs to be compared with ram only approach)            
       
        !!! H_full filling
        write(*,*) 'expected counter: =', N_spin_max*(N_spin_max+1)/2
        counter = 0
        ind_ia = 0
        do ind_i = 1, N_spin_max
            ind_ia = ind_ia + 1
            write(20,*) ind_ia
            
            do ind_j = ind_i, N_spin_max
                call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
                
                if (ind_i  == ind_j) then
                    !write(20, *)
                    write(21, *) ind_j
                    write(22, *) H_full_ij
                    counter = counter + 1
                else if (H_full_ij .NE. 0.0d0) then
                    !write(20, *) 
                    ind_ia = ind_ia + 1
                    write(21, *) ind_j
                    write(22, *) H_full_ij
                    counter = counter + 1
                else
                    counter = counter + 1
                end if
                
            end do
            
        end do
        
        write(20,*) ind_ia+1
        
        write(*,*) 'counter test = ', counter 
        no_of_nonzero = ind_ia
        write(*,*) 'number of non-zero elements = ', no_of_nonzero
        
        close(20)
        close(21)
        close(22)

    end subroutine H_XXX_feast_vec_fill
    
    subroutine H_XXX_feast_vec_diag(N_spin, no_of_nonzero)
        use omp_lib
        use mkl_vsl
        use omp_mkl_test
        implicit none
                    
        integer, intent(in) :: N_spin, no_of_nonzero
        integer :: i, j, N_spin_max, info, m_eig, n, loop, m0, fpm(128), ja_temp, ia_temp

        CHARACTER*1 :: uplo
        double precision, allocatable :: values_array(:), x(:,:), e(:), res(:)
        integer, allocatable :: ia(:), ja(:)
        double precision :: emin, emax, epsout, norm, val_arr_temp

        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        !Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

        file_name1 = 'H_eigenvals_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        !file_name2 = 'H_eigenvectors_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name3 = 'H_norm_test_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'

        open (unit=30, file= trim(file_name1), recl=512)
        !open (unit=31, file=trim(file_name2), recl=512)
        open (unit=32, file= trim(file_name3) , recl=512)

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        N_spin_max = 2**N_spin
        allocate( values_array(no_of_nonzero), ia(N_spin_max+1), ja(no_of_nonzero) )
        
        !reading of exsisting files
        open (unit=20, file="ia.dat", recl=512)
        open (unit=21, file="ja.dat", recl=512)
        open (unit=22, file="values_array.dat", recl=512)
        
        do i=1, N_spin_max+1           
            read(20, *) ia_temp
            ia(i) = ia_temp
        end do
            
        do i=1, no_of_nonzero
            read(22,*) val_arr_temp
            values_array(i) = val_arr_temp
            read(21, *) ja_temp
            ja(i) = ja_temp
        end do
        
        close(20)
        close(21)
        close(22)

        !n=n     ! Sets the size of the problem
        !a=non_zero_array     ! Array containing the nonzero elements of the upper triangular part of the matrix A
        call feastinit (fpm)  ! function specifying default parameters fpm of FEAST algorithm
        fpm(1) = 1
        fpm(2) = 12 !can be more, can be less
        fpm(3) = 8 !eps 10^-fpm(3)
        fpm(4) = 20 ! max number of feast loops
        fpm(5) = 0 !initial subspace
        fpm(6) = 0! stopping criterion
        fpm(7) = 5 !Error trace sigle prec stop crit
        fpm(14) = 0 ! standard use of feast
        fpm(27) = 1 !check input matrices
        fpm(28) = 1 !check if B is positive definite?         
        uplo='U' ! If uplo = 'U', a stores the upper triangular parts of A.
        emin = -10.0d0 ! The lower ... &
        emax =  -3.0d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = N_spin_max !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
        !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = 10000
        allocate( x(n,m0), e(m0), res(m0) )
        write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_scsrev... '
        call dfeast_scsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m_eig, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) ' dfeast_scsrev info=', info
        write(*,*) 'After  dfeast_scsrev... '
            
        if (info /= 0) then
        write(*,*) 'problem with  dfeast_scsrev, info=', info
        end if 
        
        write(*,*) ' dfeast_scsrev eigenvalues found= ', m_eig
        write(30,*) 'i-th , eigenvalue'
        do i = 1 , m_eig
            write(30,*) i, ',' , e(i)
        end do
        
        !write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        !do i=1, m_eig
           ! do j=1, n
            !    write(31,*) x(j,i)
            !end do
            !write(31,*) " "
       ! end do

        write(*,*) ' dfeast_scsrev eigenvec norm:'
        do i=1, m_eig
            norm = 0.0d0
            do j=1, n
                norm = norm + x(j,i)*(x(j,i))
            end do
            write(32,*) i, norm
        end do
                    
        deallocate(values_array, ia, ja, x, e, res)
        
        close(30)
        !close(31)
        close(32)

    end subroutine H_XXX_feast_vec_diag

end module spin_systems

program spin_code
    use spin_systems
    implicit none

    integer :: N_spin, no_of_nonzero
    double precision :: J_spin 
    character(len = 12) :: N_spin_char, J_spin_char

    ! do zrobienia
    ! zapis eigenvectora do pliku 

    ! Do Macka: 
    ! Pivoting Macierzy H w csr3? - chyba nie lepiej zrobic to znajac algorytm fizyczny
    
    ! TO DO LATER
    ! Zacznij ciac odpowiednie wektory do S_z , pytanie czy to w fortranie czy nie? 
    ! chyba w fotranie + diagonalizacja za pomocą feasta blockowych macierzy 
    ! w pythonie robisz to tnac pełną macierz H
    ! pytanie do Pawła: jak brać wycinek macierzy w H blockowym do feasta (numer eigenvalues i okno energetyczne)
    ! net cdf - binarna baza danych
   
    If(command_argument_count().NE.2) Then
        write(*,*)'Error, Only N (integer) and J (double precision) is required, program stopped'
        stop 
    endif 

    call get_command_argument(1, N_spin_char)
    call get_command_argument(2, J_spin_char)

    read(N_spin_char, *) N_spin
    read(J_spin_char, *) J_spin

    !N_spin = 4
    !J_spin = 1 

    write(*,*) 'Calculation of Heisenberg chain for N =', N_spin , 'and J = ', J_spin 

    ! diagonalization using LAPACK from intel 
    !call H_XXX_diag(N_spin, J_spin)

    ! diagonalization via FEAST algorithm (saving CSR matrices to file right now)
    call sparse_dfeast_test()
    call H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero)
    call H_XXX_feast_vec_diag(N_spin, no_of_nonzero)

end program spin_code