program readdata
implicit none 


integer :: line_count1 , ios, j, i, it, resultindex, msize1, msize2,n, flops,nx,ny
integer,  dimension(:), allocatable :: X1_cord, Y1_cord, X_cord_new, Y_cord_new, uniqueuncomp, unique
real(kind=8), dimension(:), allocatable :: valxy1, valxy_new,yvec,xvec
real(kind=8) :: r
character(256) :: filename1
filename1 = "MATRIX_125_25p"
line_count1 = 0
resultindex =1
msize1= 128
flops=0
!reads the amount of lines in the input of the first matrix 
open(1, file=filename1 , status="OLD")
do 
    read(1, *, iostat=ios)
    if (ios/=0) exit
    line_count1 = line_count1 + 1 
end do
close(1)

!allocate memory and reads matrix into 3 1D arrays 
allocate(X1_cord(line_count1))
allocate(Y1_cord(line_count1))
allocate(valxy1(line_count1))

open(1, file=filename1, status="OLD")
do i = 1,line_count1
    read(1,*) Y1_cord(i), X1_cord(i), valxy1(i)
end do
close(1)


n = msize1 * (msize1 + 1) / 2
allocate(yvec(msize1)) !temporary arrays that will be used for each row/column pair
allocate(xvec(msize1))
allocate(X_cord_new(n))
allocate(Y_cord_new(n))
allocate(valxy_new(n))
allocate(uniqueuncomp(line_count1))
uniqueuncomp=0
j = 0
n = 1
do i=1, line_count1   !find all the unique indexes
   nx =0
   ny =0
   do it=1, line_count1
      if (uniqueuncomp(it) == Y1_cord(i)) then
          ny = 1 
      end if    
      if (uniqueuncomp(it) == X1_cord(i)) then
          nx = 1 
      end if
      if ((nx ==1) .AND. (ny==1)) then
              exit
      end if
    end do

    if (uniqueuncomp(it) ==0) then
       if (Y1_cord(i) == (X1_cord(i))) then
               uniqueuncomp(Y1_cord(i)) = Y1_cord(i)
               j = j +1 
               flops = flops + 1
               cycle
       end if
       if (ny == 0) then
           uniqueuncomp(Y1_cord(i)) = Y1_cord(i)
           j = j +1
           flops = flops + 1
        end if
        if (nx == 0) then 
              uniqueuncomp(X1_cord(i)) = X1_cord(i)
              j = j +1 
              flops = flops + 1
       
        end if
    end if 
end do
allocate(unique(j))

n = j
j = 1
do i=1, line_count1
   if (uniqueuncomp(i) == 0) then
           cycle
   end if
   unique(j) = uniqueuncomp(i)
   j = j + 1
   flops = flops + 1 
end do

do i=1, n
   yvec = 0
   do it=1, line_count1
      if (Y1_cord(it) == unique(i)) then
              yvec(X1_cord(it)) = valxy1(it)
      end if
      if (X1_cord(it) == unique(i)) then
              yvec(Y1_cord(it)) = valxy1(it)
      end if
   end do

   do j=1, n
      if (j < i) then
          cycle
      end if
      xvec = 0
      do it=1, line_count1
         if (Y1_cord(it) == unique(j)) then
                 xvec(X1_cord(it)) = valxy1(it)
         end if
         if (X1_cord(it) == unique(j)) then
                 xvec(Y1_cord(it)) = valxy1(it)
         end if
      end do
      r = 0.0
      do it=1, msize1
         r = r + yvec(it) * xvec(it)
         flops = flops + 1
      end do 
      if (r == 0.0) then
        cycle
      end if
      
      valxy_new(resultindex) = r
      X_cord_new(resultindex) = j
      Y_cord_new(resultindex) = i
      resultindex = resultindex +1 
      flops = flops + 1
   end do
end do


print *,valxy_new
print *,Y_cord_new
print *,X_cord_new
print*,flops

end program readdata
