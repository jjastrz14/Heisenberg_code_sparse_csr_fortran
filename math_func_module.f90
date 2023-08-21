module math_functions
    implicit none
    contains

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

end module math_functions