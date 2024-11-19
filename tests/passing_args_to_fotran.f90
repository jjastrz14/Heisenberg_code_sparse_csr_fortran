! example.f90
program main
    implicit none
    character(len=*), parameter :: VERSION = '1.0'
    character(len=32)           :: arg
    integer                     :: i

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)

        select case (arg)
            case ('-v', '--version')
                print '(2a)', 'version ', VERSION
                stop

            case ('-h', '--help')
                call print_help()
                stop

            case default
                print '(2a, /)', 'unrecognised command-line option: ', arg
                call print_help()
                stop
        end select
    end do

    print '(a)', 'Hello, World!'
contains
    subroutine print_help()
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -v, --version     print version information and exit'
        print '(a, /)', '  -h, --help        print usage information and exit'
    end subroutine print_help
end program main

! compile as 
! $ gfortran13 -o example example.f90
! $ ./example --help
! command-line options:

 ! -v, --version     print version information and exit
 ! -h, --help        print usage information and exit

!$ ./example --version
!version 1.0