program multiply_matrices
    use dense_module
    use sparse_module
    use shared_module

    implicit none

    ! Declare variables and matrices
    integer :: size, Nmult, nnz1, nnz2
    real :: totalTime
    character(len=100) :: fileA, fileB
    real(kind=8), allocatable :: A(:,:), B(:,:), C(:,:), Val1(:), Val2(:)
    integer, allocatable :: Row1(:), Col1(:), Row2(:), Col2(:)
    character(len=10) :: size_str, debug_flag_str
    logical :: debug_flag

    ! Take arguments
    if (command_argument_count() < 3) then
        print *, "Usage: dense_mult.f90 fileA fileB size [debug_flag]"
        stop
    end if
    call get_command_argument(1, fileA)
    call get_command_argument(2, fileB)
    call get_command_argument(3, size_str)
    read(size_str, *) size

    ! Optional debugging flag
    debug_flag = .false.
    if (command_argument_count() >= 4) then
        call get_command_argument(4, debug_flag_str)
        if (debug_flag_str /= "0") debug_flag = .true.
    end if

    ! Allocate matrices and set to 0
    allocate(A(size, size), B(size, size), C(size, size))

    ! Read matrices from files and store in 2D dense arrays
    call store_dense(A, size, fileA)
    call store_dense(B, size, fileB)

    ! Generate sparse matrices
    call get_nnz(A, size, nnz1)
    allocate(Row1(nnz1), Col1(nnz1), Val1(nnz1))
    call get_nnz(B, size, nnz2)
    allocate(Row2(nnz2), Col2(nnz2), Val2(nnz2))
    call dense2sym_sparse(Row1, Col1, Val1, A, size)
    call dense2sym_sparse(Row2, Col2, Val2, B, size)

    ! Debugging: Print matrices if debug_flag is true
    if (debug_flag) then
        print *, "Matrix A:"
        call print_matrix(size, A)
        print *
        print *, "Matrix B:"
        call print_matrix(size, B)
        print *
        print *, "#################################"
    end if

    ! Measure and print multiplication times for each dense method
    call measure_mult("dense naive", A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2, debug_flag)

    call measure_mult("sparse symmetry", A, B, C, size, Nmult, totalTime, &
    Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2, debug_flag)

    print *, "Warning! # Multiplications will be 0 as cannot directly measure LAPACK's DGEMM operations"
    call measure_mult("LAPACK", A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2, debug_flag)

    deallocate(A, B, C, Row1, Col1, Val1, Row2, Col2, Val2)
end program multiply_matrices
