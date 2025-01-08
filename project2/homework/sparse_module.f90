Module sparse_module
    implicit none
    Contains

    function get_nlines(file, size) result(nlines)
        implicit none
        character(len=*), intent(in) :: file
        integer, intent(in) :: size
        integer :: nlines

        integer :: io, ios
        real(kind=8) :: value
        integer :: row_counts(size)

        open(unit=io, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Cannot open file ", trim(file)
            stop
        end if

        nlines = 0
        do
            read(io, *, iostat=ios)
            nlines = nlines + 1 
            if (ios < 0) exit       ! End of file
            if (ios > 0) then       ! Error reading the line
                print *, "Error: Problem reading the file!"
                stop
            end if
        end do
        close(io)
    end function get_nlines

    subroutine store_sparse(size, file, values, row_ind, col_ind) ! Store matrix explicitly in 2D array
        implicit none
        integer, intent(in) :: size
        character(len=*), intent(in) :: file
        real(kind=8), intent(inout) :: values(:)
        integer, intent(inout) :: row_ind(:), col_ind(:)

        integer :: i, j, n, io, ios          
        real(kind=8) :: value    

        ! Checks if file exists, otherwise it stops the program with error
        open(unit=io, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Cannot open file ", trim(file)
            stop
        end if

        row_ind = 0
        col_ind = 0
        values = 0.0d0
        n = 1
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
                row_ind(n) = i
                col_ind(n) = j
                values(n) = value
                n = n + 1
            else
                print *, "Warning: For file", trim(file),", indices out of bounds for size=", size,": i=", i, ", j=", j
                print *, "Ensure both matrices are of size <=", size
                stop
            end if
        end do
    end subroutine store_sparse

    subroutine mult_sparse_sim(Row1, Col1, Val1, Row2, Col2, Val2, C, size, nnz1, nnz2, Nmult) !Multiplies two sparse matrices in the symmtrized forms and outputs a dense matrix
        implicit none
        integer, intent(in):: size, nnz1, nnz2
        integer, intent(inout):: Nmult
        integer, intent(in) :: Row1(:), Col1(:), Row2(:), Col2(:) 
        real(kind=8), intent(in) :: Val1(:), Val2(:)
        real(kind=8), intent(inout) :: C(:,:) 

        integer :: i, j

        C = 0.d0

        ! Thanks Jacopo
        do i = 1,nnz1
            do j = 1,nnz2
                if (Col1(i).eq.Row2(j)) then
                    C(int(Row1(i)),int(Col2(j))) = C(int(Row1(i)),int(Col2(j))) + Val1(i)*Val2(j)
                    Nmult = Nmult + 1
                elseif (Col1(i).eq.Col2(j).and. Col2(j).ne.Row2(j)) then
                    C(int(Row1(i)),int(Row2(j))) = C(int(Row1(i)),int(Row2(j))) + Val1(i)*Val2(j)
                    Nmult = Nmult + 1
                endif
            enddo
        enddo
        do i = 1,nnz1
            do j = 1,nnz2
                if (Row1(i).eq.Col2(j).and. Row1(i).ne.Col1(i).and.Col2(j).ne.Row2(j)) then
                    C(int(Col1(i)),int(Row2(j))) = C(int(Col1(i)),int(Row2(j))) + Val1(i)*Val2(j)
                    Nmult = Nmult + 1
                elseif (Row1(i).eq.Row2(j).and. Row1(i).ne.Col1(i)) then
                    C(int(Col1(i)),int(Col2(j))) = C(int(Col1(i)),int(Col2(j))) + Val1(i)*Val2(j)
                    Nmult = Nmult + 1
                endif
            enddo
        enddo

    end subroutine mult_sparse_sim

    subroutine sym_sparse2dense(Row, Col, Val, A, nnz)
        implicit none
        integer, intent(in):: nnz
        integer, intent(in) :: Row(:), Col(:)
        real(kind=8), intent(in) :: Val(:)
        real(kind=8), intent(inout) :: A(:,:) 

        integer :: i, j

        A = 0.d0

        do i = 1, nnz
            do j = 1, nnz
                A(Row(i),Col(j)) = Val(j)
                A(Col(i),Row(j)) = Val(j) ! No need for if in diagonal as it will override the same value
            end do
        end do
    end subroutine sym_sparse2dense

    subroutine dense2sym_sparse(Row, Col, Val, A, size)
        implicit none
        integer, intent(in):: size
        real(kind=8), intent(in) :: A(:,:) 
        integer, intent(inout) :: Row(:), Col(:)
        real(kind=8), intent(inout) :: Val(:)

        integer :: i, j, n
        
        Row = 0
        Col = 0
        Val = 0.d0
        n = 1

        do j = 1, size
            do i = j, size
                if (A(i,j).ne.0) then
                    Row(n) = i
                    Col(n) = j
                    Val(n) = A(i,j)
                    n = n + 1
                end if
            end do
        end do
    end subroutine dense2sym_sparse

    subroutine get_nnz(A, size, nnz)
        implicit none
        integer, intent(in):: size
        real(kind=8), intent(in) :: A(:,:) 
        integer, intent(inout) :: nnz

        integer :: i, j
        
        nnz = 0

        do j = 1, size
            do i = j, size
                if (A(i,j).ne.0) then
                    nnz = nnz + 1
                end if
            end do
        end do
    end subroutine get_nnz

end Module sparse_module

!    subroutine sort_3arr_by1(array1, array2, array3, sorted1, sorted2, sorted3)
!        implicit none !Sorts 3 arrays by the first given array
!        real, intent(in) :: array1(:), array2(:), array3(:)
!        real, intent(inout) :: sorted1(:), sorted2(:), sorted3(:)
!        integer :: i, j, n
!        real :: temp
!
!        ! Copy input arrays to the sorted arrays
!        sorted1 = array1
!        sorted2 = array2
!        sorted3 = array3
!
!        ! Bubble sort based on sorted1
!        do i = 1, n - 1
!            do j = 1, n - i
!                if (sorted1(j) > sorted1(j + 1)) then
!                    ! Swap sorted1
!                    temp = sorted1(j)
!                    sorted1(j) = sorted1(j + 1)
!                    sorted1(j + 1) = temp
!
!                    ! Swap sorted2
!                    temp = sorted2(j)
!                    sorted2(j) = sorted2(j + 1)
!                    sorted2(j + 1) = temp
!
!                    ! Swap sorted3
!                    temp = sorted3(j)
!                    sorted3(j) = sorted3(j + 1)
!                    sorted3(j + 1) = temp
!                end if
!            end do
!        end do
!    end subroutine sort_three_arrays
!
!    function count_repeats(array, idx, val) result(count)
!            implicit none
!            integer, intent(in) :: array(:) 
!            integer, intent(in) :: idx      
!            integer, intent(in) :: val      
!            integer :: count                
!            integer :: i
!
!            count = 0
!            ! Count occurrences of val starting from the given index
!            i = idx
!            do 
!                if (array(i) == val) then
!                    count = count + 1
!                    i = i + 1
!                else
!                    exit
!                end if
!            end do
!    end function count_repeats

