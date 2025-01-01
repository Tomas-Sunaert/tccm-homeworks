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

!Debugging
print *, 'Matrix A:'
call print_matrix(size,A)
print *
print *, 'matrix B:'
call print_matrix(size,B)
print *
print *, "#################################"

    ! Measure and print multiplication times for each dense method
    call measure_mult("dense naive", A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2)
print *, 'matrix C:' ! Debugging
call print_matrix(size,C)
print *
print *, "#################################"

    call measure_mult("sparse symmetry", A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2)
print *, 'matrix C:' !Debugging
call print_matrix(size,C)
print *
print *, "#################################"

print *, "Warning! # Multiplications will be 0 as cannot directly measure LAPACK\'s DGEMM operations"
    call measure_mult("LAPACK", A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2)
print *, 'matrix C:' !Debugging
call print_matrix(size,C)
print *

    deallocate(A, B, C, Row1, Col1, Val1, Row2, Col2, Val2)
end program multiply_matrices