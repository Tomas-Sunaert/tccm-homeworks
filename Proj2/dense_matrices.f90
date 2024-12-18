program multiply_matrices
    use dense_module
    implicit none

    ! Declare variables and matrices
    integer :: size, Nmult
    real :: totalTime
    character(len=100) :: fileA, fileB
    real(kind=8), allocatable :: A(:,:), B(:,:), C(:,:)
    character(len=10) :: size_str

    !Take arguments  
    if (command_argument_count() /= 3) then
        print *, "Usage: dense_mult.f90 fileA fileB size"
        stop
    end if
    call get_command_argument(1, fileA)
    call get_command_argument(2, fileB)
    call get_command_argument(3, size_str)
    read(size_str, *) size
    ! Allocate matrices
    allocate(A(size, size), B(size, size), C(size, size))

    ! Read matrices from files and store in 2D dense arrays
    call store_dense(A, size, fileA)
    call store_dense(B, size, fileB)

    ! Measure and print multiplication times for each dense method
    call measure_mult("naive", A, B, C, size, Nmult, totalTime)
    call measure_mult("symmetry", A, B, C, size, Nmult, totalTime)
    call measure_mult("LAPACK", A, B, C, size, Nmult, totalTime)

    deallocate(A,B,C)
end program multiply_matrices