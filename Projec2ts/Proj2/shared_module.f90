Module shared_module
    use dense_module
    use sparse_module
    use dgemm_module
    implicit none
    Contains

    subroutine debug_print(C, size, method, debug_flag)
        logical, intent(in)  :: debug_flag
        integer, intent(in) :: size
        character(len=*), intent(in) :: method      
        real(kind=8), intent(in) :: C(:,:)

        if (debug_flag) then
            print *, "Matrix C after", method,":"
            call print_matrix(size, C)
            print *
            print *, "#################################"
        end if
    end subroutine debug_print

    subroutine call_mult(method, A, B, C, size, Nmult, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2)
        implicit none
        character(len=*), intent(in) :: method
        real(kind=8), intent(in) :: A(:,:), B(:,:), Val1(:), Val2(:)
        real(kind=8), intent(out) :: C(:,:)
        integer, intent(in) :: size, Row1(:), Row2(:), Col1(:), Col2(:), nnz1, nnz2
        integer, intent(out) :: Nmult
        
        ! Call multiplication method
        if (method == "dense naive") then
            call mult_dense_naive(A, B, C, size, Nmult)
        elseif (method == "LAPACK") then
            C = 0.0d0
            call dgemm('N', 'N', size, size, size, 1.0d0, A, size, B, size, 0.0d0, C, size,Nmult)
        elseif (method == "sparse symmetry") then
            call mult_sparse_sim(Row1, Col1, Val1, Row2, Col2, Val2, C, size, nnz1, nnz2, Nmult)
        else
            print *, "Error: Invalid method specified: ", method
            stop
        end if
    end subroutine call_mult
      
    subroutine measure_mult(method, A, B, C, size, Nmult, totalTime, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2, debug_flag) ! Subroutine to measure matrix multiplication time
        implicit none
        character(len=*), intent(in) :: method
        real(kind=8), intent(in) :: A(:,:), B(:,:), Val1(:), Val2(:)
        real(kind=8), intent(out) :: C(:,:)
        integer, intent(in) :: size, Row1(:), Row2(:), Col1(:), Col2(:), nnz1, nnz2
        integer, intent(out) :: Nmult
        real, intent(out) :: totalTime
        logical, intent(in)  :: debug_flag

        integer :: t
        real(kind=8) :: beginning, end

        print *, "Matrix multiplication using ", method, " method:"
        Nmult = 0
        totalTime = 0

        call cpu_time(beginning)
        call call_mult(method, A, B, C, size, Nmult, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2) ! Call multiplication method
        call cpu_time(end)
        totalTime = real(end - beginning)

        ! In case one matrix multiplication was too slow to be measured
        if (totalTime == 0) then
            call cpu_time(beginning)
            Nmult = 0
            do t = 1, 1000
                call call_mult(method, A, B, C, size, Nmult, Row1, Row2, Col1, Col2, Val1, Val2, nnz1, nnz2)
            end do
            call cpu_time(end)
            Nmult = Nmult / 1000
            totalTime = real(end - beginning) / 1000
        end if

        print *, "Wall time: ", totalTime, "s"
        print *, "# Multiplications: ", Nmult
        print *
        call debug_print(C, size, method, debug_flag)
    end subroutine measure_mult
end Module shared_module
