module matrix_io
  implicit none
contains
  subroutine read_mat(file_path, row, col, val, num_lines, dimensions, dense)
    implicit none

    character(len=*), intent(in) :: file_path
    real*8, allocatable, intent(out) :: row(:), col(:), val(:)
    real*8, allocatable, intent(out) :: dense(:,:)
    integer, intent(out) :: num_lines, dimensions

    integer :: i, io_status, max_row, max_col, max_index
    real*8 :: row_temp, col_temp, val_temp

    open(unit=13, file=file_path, status="old", action="read")
    ! Conta il numero di linee
    num_lines = 0
    do
      read(13, *, iostat=io_status)
      if (io_status /= 0) exit
      num_lines = num_lines + 1
    end do

    rewind(13)

    ! Alloca gli array
    allocate(row(num_lines), col(num_lines), val(num_lines))
    row = 0.d0
    col = 0.d0
    val = 0.d0

    ! Legge i valori
    do i = 1, num_lines
      read(13, *) row(i), col(i), val(i)
    end do
    close(13)

    ! Trova il massimo indice
    max_row = maxval(row)
    max_col = maxval(col)
    max_index = max(max_row, max_col)

    ! Alloca la matrice densa
    allocate(dense(max_index, max_index))
    dense = 0.d0
    do i = 1, num_lines
      dense(int(row(i)), int(col(i))) = val(i)
      dense(int(col(i)), int(row(i))) = val(i) ! Assumendo matrice simmetrica
    end do


    dimensions = max_index
  end subroutine read_mat
  subroutine SMaMul(Row1,Row2,Col1,Col2,Val1,Val2, Res, size1, size2, dimension1, dimension2)
  implicit none
  integer :: dimension1, dimension2, size1, size2, big_dim
  integer :: i,j, n_operations=0
  real*8,allocatable :: Row1(:), Col1(:), Val1(:) ! FIRST MATRIX DESCRIPPTORS
  real*8,allocatable :: Row2(:), Col2(:), Val2(:) ! SECOND MATRIX DESCRIPTORS
  real*8, allocatable :: Res(:,:) ! RESULT MATRIX
  real*8 :: Mel ! matrix element (for accumulation)

  write(*,*) 'Dim 1=', dimension1
  write(*,*) 'Dim 2=', dimension2
  
  write(*,*) 'Size 1=', size1
  write(*,*) 'Size 2=', size2
  ! check that we are trying to multuply to matrices that have same dimentionality
  if (dimension1.ne.dimension2) then
   write(*,*) 'Error, Matrix 1 and Matrix 2 have different dimentions!'
   stop
  endif

  !allocate(Row1(dimension1),Col1(dimension1),Val1(dimension1))
  !allocate(Row1(dimension2),Col1(dimension2),Val1(dimension2))


 
  big_dim=max(dimension1,dimension2)
  allocate(Res(big_dim,big_dim))



  Res = 0.d0 !set every element to 0

  do i = 1,size1
      do j = 1,size2
          if (Col1(i).eq.Row2(j)) then
              Res(int(Row1(i)),int(Col2(j))) = Res(int(Row1(i)),int(Col2(j))) + Val1(i)*Val2(j)
          elseif (Col1(i).eq.Col2(j).and. Col2(j).ne.Row2(j)) then
              Res(int(Row1(i)),int(Row2(j))) = Res(int(Row1(i)),int(Row2(j))) + Val1(i)*Val2(j)
          endif
      enddo
  enddo
write(*,*) 'N 2'

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

integer :: size1, size2, dimensions1, dimensions2
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


call read_mat(file_path1, row1, col1, val1, size1, dimensions1, dense1)
call read_mat(file_path2, row2, col2, val2, size2, dimensions2, dense2)

allocate(ResDense(dimensions1,dimensions1))
ResDense=matmul(dense1,dense2)
call SMaMul(Row1,Row2,Col1,Col2,Val1,Val2, Res, size1, size2, dimensions1, dimensions2)







open(19,file='MaMul')
do i = 1,dimensions1
    write(19,*) ResDense(i,:)
enddo
close(19)

open(19,file='My_MaMul')
do i = 1,dimensions1
    write(19,*) Res(i,:)
enddo






end program

