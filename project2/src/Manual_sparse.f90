! this is this is the easy way: where we 

program main
implicit none

real*8,allocatable :: Row(:), Col(:), Val(:) !3 allocatable arrays 
real*8, allocatable :: Dense(:,:) ! the matrix representation

integer :: i,j,k,l ! indexes
integer :: num_lines ! number of lines
integer :: io_status ! error checking during the execution
integer :: max_row, max_col, max_index

character*64 :: file_path, file_name

! ask the name of the file
write(*,*) 'Choose:'
CALL SYSTEM('ls ../data')
write(*,*) 'Name of the file with the data?'
read(*,*) file_name

! setting the path right
file_path = '../data/'//file_name


open(unit=13, file=file_path, status="old", action="read")
! counting the number of lines
num_lines = 0
do
 read(13, *, iostat=io_status)
  if (io_status /= 0) exit
   num_lines = num_lines + 1
end do

REWIND(13)

! now i can allocate the arrays 
allocate(Row(num_lines), Col(num_lines), Val(num_lines))

! inizialization for good practice
Row = 0.d0
Col = 0.d0
Val = 0.d0


! reading the value
do i = 1, num_lines
 read(13,*) Row(i), Col(i), Val(i)
enddo
close(13)

max_row = maxval(Row)
max_col = maxval(Col)


max_index = max(max_row,max_col)

! now that we know what is the max index that occurs, we just allocate a big enough dense 2D matrix

allocate(Dense(max_index,max_index))

Dense=0.d0
do i = 1, num_lines
 Dense(int(Row(i)),int(Col(i))) = Val(i)
 Dense(int(Col(i)),int(Row(i))) = Val(i)
enddo

open(12, file='Dense_matrix')


do i = 1, max_index
  write(12,*) Dense(i,:)
enddo

close(12)




call system('gnuplot -e "set xlabel ''Column index''; set ylabel ''Row index'';  plot ''Dense_matrix'' matrix with image; pause -1"')







end program
