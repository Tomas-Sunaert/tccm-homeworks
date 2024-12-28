Module dense_module
    implicit none
    Contains

    subroutine print_matrix(size,A) ! Print matrix to screen
        implicit none
        integer, intent(in) :: size
        real(kind=8), intent(in) :: A(size, size)

        integer :: i

        do i = 1, size
            print *, A(i, 1:size)
        end do
    end subroutine print_matrix

    subroutine store_dense(matrix, size, file) ! Store matrix explicitly in 2D array
        implicit none
        real(kind=8), intent(inout) :: matrix(:,:)
        integer, intent(in) :: size
        character(len=*), intent(in) :: file

        integer :: i, j, io, ios          
        real(kind=8) :: value    

        ! Checks if file exists, otherwise it stops the program with error
        open(unit=io, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Cannot open file ", trim(file)
            stop
        end if

        matrix = 0.0d0
        ! Read file line by line and store to simetrized matrix
        do
            read(io, *, iostat=ios) i, j, value
            if (ios < 0) then
                exit                   ! End of file
            
            else if (ios > 0) then                   ! Error reading the line
                print *, "Error: Problem reading the file. Format must be: i(int) j(int) value(double)"
                stop
            
            ! Check indices are within bounds
            else if (i >= 1 .and. i <= size .and. j >= 1 .and. j <= size) then
                matrix(i, j) = value
                matrix(j, i) = value !Matrices are simmetric
            else
                print *, "Warning: For file", trim(file),", indices out of bounds for size=", size,": i=", i, ", j=", j
                print *, "Ensure both matrices are of size <=", size
                stop
            end if
        end do

        close(io)
        
    end subroutine store_dense

    subroutine mult_dense_naive(A, B, C, size, Nmult) !Naive matrix multiplication
        implicit none
        real(kind=8), intent(in) :: A(:,:), B(:,:)
        real(kind=8), intent(inout) :: C(:,:)
        integer, intent(in) :: size        
        integer, intent(inout) :: Nmult

        integer :: i, j, k

        C = 0.0d0
        do i = 1, size
            do j = 1, size
                do k = 1, size
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                    Nmult = Nmult + 1
                end do
            end do
        end do        
    end subroutine mult_dense_naive

    subroutine mult_dense_sim(A, B, C, size, Nmult) ! Use symetries to perform les operations
        implicit none
        real(kind=8), intent(in) :: A(:,:), B(:,:)
        real(kind=8), intent(inout) :: C(:,:)
        integer, intent(in) :: size        
        integer, intent(inout) :: Nmult

        integer :: i, j, k

        C = 0.0d0
        do i = 1, size
            do j = i, size ! Only iterate over the upper triangular part (j >= i) 
                do k = 1, size
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                    Nmult = Nmult + 1
                end do
                C(j, i) = C(i, j) !This is not true
            end do
        end do
    end subroutine mult_dense_sim

    subroutine call_mult(method, A, B, C, size, Nmult)
        implicit none
        character(len=*), intent(in) :: method
        real(kind=8), intent(in) :: A(:,:), B(:,:)
        real(kind=8), intent(out) :: C(:,:)
        integer, intent(in) :: size
        integer, intent(out) :: Nmult

        ! Call multiplication method
        if (method == "naive") then
            call mult_dense_naive(A, B, C, size, Nmult)
        elseif (method == "symmetry") then
            call mult_dense_sim(A, B, C, size, Nmult)
        elseif (method == "LAPACK") then
            C = 0.0d0
            call dgemm('N', 'N', size, size, size, 1.0d0, A, size, B, size, 0.0d0, C, size)
        else
            print *, "Error: Invalid method specified: ", method
            stop
        end if
    end subroutine call_mult
      
    subroutine measure_mult(method, A, B, C, size, Nmult, totalTime) ! Subroutine to measure matrix multiplication time
        implicit none
        character(len=*), intent(in) :: method
        real(kind=8), intent(in) :: A(:,:), B(:,:)
        real(kind=8), intent(out) :: C(:,:)
        integer, intent(in) :: size
        integer, intent(out) :: Nmult
        real, intent(out) :: totalTime
        integer :: beginning, end, rate, t

        print *, "Matrix multiplication using ", method, " method:"
        Nmult = 0
        call system_clock(beginning, rate)
        call call_mult(method, A, B, C, size, Nmult) ! Call multiplication method
        call system_clock(end)
        totalTime = real(end - beginning) / real(rate)

        ! In case one matrix multiplication was too slow to be measured
        if (totalTime == 0) then
            call system_clock(beginning, rate)
            Nmult = 0
            do t = 1, 1000
                call call_mult(method, A, B, C, size, Nmult)
            end do
            call system_clock(end)
            Nmult = Nmult / 1000
            totalTime = real(end - beginning) / real(rate) / 1000
        end if

        print *, "Wall time: ", totalTime, "s"
        print *, "# Multiplications: ", Nmult
        print *
    end subroutine measure_mult

end Module dense_module

