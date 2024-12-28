! here is a module containing all the subroutines and global variable (none) in this script
module matrix_io
  implicit none
contains
  subroutine read_mat(file_path, row, col, val, num_lines, dimensions, dense)
    implicit none

    character(len=*), intent(in) :: file_path !path
    real*8, allocatable, intent(out) :: row(:), col(:), val(:) !arrays in wich to store the 3 1D arrays
    real*8, allocatable, intent(out) :: dense(:,:) !2D dense format
    integer, intent(out) :: num_lines, dimensions ! array elements of the 1D arrays as well as the max index 

    integer :: i, io_status, max_row, max_col, max_index !error check and utilities
    real*8 :: row_temp, col_temp, val_temp ! utilities

    open(unit=13, file=file_path, status="old", action="read")
    ! Counts the lines
    num_lines = 0
    do
      read(13, *, iostat=io_status)
      if (io_status /= 0) exit
      num_lines = num_lines + 1
    end do

    rewind(13)

    ! Now that i know how many lines it takes i allocate the arrays
    allocate(row(num_lines), col(num_lines), val(num_lines))
    ! Initialization (safety)
    row = 0.d0
    col = 0.d0
    val = 0.d0

    ! Reading the values from the input file
    do i = 1, num_lines
      read(13, *) row(i), col(i), val(i)
    end do
    close(13)

    ! Finds the max index value in the matrix
    max_row = maxval(row)
    max_col = maxval(col)
    max_index = max(max_row, max_col)

    ! Allocates a 2D array for the dense rappresentation
    allocate(dense(max_index, max_index))
    dense = 0.d0
    do i = 1, num_lines
      dense(int(row(i)), int(col(i))) = val(i)
      dense(int(col(i)), int(row(i))) = val(i) ! Taking care of the symetry
    end do


    dimensions = max_index ! For the output, i assign to dimentions the value of max_index


  end subroutine read_mat



  ! NOW I FOCUS ON THE REAL DEAL
  subroutine SMaMul(Row1,Row2,Col1,Col2,Val1,Val2, Res, size1, size2, dimension1, dimension2)
  implicit none
  integer :: dimension1, dimension2, size1, size2, big_dim
  integer :: i,j, n_operations=0
  real*8,allocatable :: Row1(:), Col1(:), Val1(:) ! FIRST MATRIX DESCRIPPTORS
  real*8,allocatable :: Row2(:), Col2(:), Val2(:) ! SECOND MATRIX DESCRIPTORS
  real*8, allocatable :: Res(:,:) ! RESULT MATRIX
  real*8 :: Mel ! matrix element (for accumulation)
  
  ! Shows the matrix dimentions
  write(*,*) 'Dim 1=', dimension1
  write(*,*) 'Dim 2=', dimension2
  ! And relative sizes of the 1D arrays we are working with (both for debugging and because is interesting to know)
  write(*,*) 'Size 1=', size1
  write(*,*) 'Size 2=', size2
  ! check that we are trying to multuply to matrices that have same dimentionality
  if (dimension1.ne.dimension2) then
   write(*,*) 'Error, Matrix 1 and Matrix 2 have different dimentions!'
   stop
  endif

  !allocate(Row1(dimension1),Col1(dimension1),Val1(dimension1))
  !allocate(Row1(dimension2),Col1(dimension2),Val1(dimension2))


  ! find the minimum size to allocate the mamori with big_dim
  big_dim=max(dimension1,dimension2)
  allocate(Res(big_dim,big_dim)) ! and then allocate the 2D form



  Res = 0.d0 !set every element to 0, in order to prepare for the accmulation
  
  ! I prototyped this thing on MathLab and got it right, try to understand it
  ! All this is just not to rely on new arrays that would waste memory
  do i = 1,size1
      do j = 1,size2
          if (Col1(i).eq.Row2(j)) then
              Res(int(Row1(i)),int(Col2(j))) = Res(int(Row1(i)),int(Col2(j))) + Val1(i)*Val2(j)
          elseif (Col1(i).eq.Col2(j).and. Col2(j).ne.Row2(j)) then
              Res(int(Row1(i)),int(Row2(j))) = Res(int(Row1(i)),int(Row2(j))) + Val1(i)*Val2(j)
          endif
      enddo
  enddo

  do i = 1,size1
      do j = 1,size2
          if (Row1(i).eq.Col2(j).and. Row1(i).ne.Col1(i).and.Col2(j).ne.Row2(j)) then
              Res(int(Col1(i)),int(Row2(j))) = Res(int(Col1(i)),int(Row2(j))) + Val1(i)*Val2(j)
          elseif (Row1(i).eq.Row2(j).and. Row1(i).ne.Col1(i)) then
              Res(int(Col1(i)),int(Col2(j))) = Res(int(Col1(i)),int(Col2(j))) + Val1(i)*Val2(j)
          endif
      enddo
  enddo

return



  end subroutine
end module matrix_io


program main
use matrix_io
implicit none

integer :: size1, size2, dimensions1, dimensions2, big_enough
integer :: i,j

real*8,allocatable :: Row1(:), Col1(:), Val1(:) ! FIRST MATRIX DESCRIPPTORS
real*8,allocatable :: Row2(:), Col2(:), Val2(:) ! SECOND MATRIX DESCRIPTORS
real*8, allocatable :: Dense1(:,:) ! FIRST MATRIX IN DENSE FORM
real*8, allocatable :: Dense2(:,:) ! SECOND MATRIX IN DENSE FORM
real*8, allocatable :: Res(:,:), ResDense(:,:)


character*64 :: file_path1, file_path2 ! THE 2 PATHS TO THE DATA FOLDER
character*64 :: file_name1, file_name2 ! NAMES FOR THE 2 MATRICES

! ask the name of the file
write(*,*) 'Choose first matrix:'
CALL SYSTEM('ls ../data')
write(*,*) 'Name of the file with the data?'
read(*,*) file_name1

! setting the path right
file_path1 = '../data/'//file_name1

! ask the name of the file
write(*,*) 'Choose second matrix:'
CALL SYSTEM('ls ../data')
write(*,*) 'Name of the file with the data?'
read(*,*) file_name2

! setting the path right
file_path2 = '../data/'//file_name2

! i call the suberoutines i made to read the matrices in sparse 1D arrays form
call read_mat(file_path1, row1, col1, val1, size1, dimensions1, dense1)
call read_mat(file_path2, row2, col2, val2, size2, dimensions2, dense2)



! I allocate a matrix with dimentions big enough
big_enough=max(dimensions1,dimensions2)
allocate(ResDense(big_enough,big_enough))
! compute the matrix porduct with the matmul function for reference
ResDense=matmul(dense1,dense2)

! compute with MY subroutine the matrix product
call SMaMul(Row1,Row2,Col1,Col2,Val1,Val2, Res, size1, size2, dimensions1, dimensions2)






! write the intrinsic function reference in MaMul
open(19,file='MaMul')
do i = 1,dimensions1
    write(19,*) ResDense(i,:)
enddo
close(19)


! write mu result in My_MaMul
open(19,file='My_MaMul')
do i = 1,dimensions1
    write(19,*) Res(i,:)
enddo






end program

