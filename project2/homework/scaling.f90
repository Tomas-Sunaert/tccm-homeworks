Program scaling
        use dense_module
        use sparse_module
        use dgemm_module
        use genmat_module
        use shared_module
        implicit none 
        
        integer :: size, sizemin, sizemax, density,densitymin,densitymax, Nmult,stepsize, nnz1, nnz2, j,i, x
        integer, allocatable ::  Row1(:), Row2(:), Col1(:),Col2(:)
        real(kind=8), allocatable :: A(:,:), B(:,:), C(:,:), Xsize(:),Nsize(:),tsize(:),Xdens(:), Ndens(:), &
        tdens(:), Val1(:), Val2(:)
        real :: totalTime
        character(len=10) :: method
        log10ical :: debug_flag
        debug_flag = .false.
        sizemin = 5
        sizemax = 500
        stepsize = 15
        densitymin = 1
        densitymax = 50
        size = sizemin
        density = 25
        method = 'sparse'
        i = 1
        x= int(real(sizemax-sizemin)/stepsize)
        allocate(Xsize(x),tsize(x),Nsize(x))
        do while(size < sizemax)
           allocate(A(size,size),B(size,size),C(size,size)) 
           Xsize(i) =  size
           Nmult = 0 
           totalTime = 0 
           call genmat(A,size,density)
           call genmat(B,size,density)
           call get_nnz(A,size,nnz1)
           call get_nnz(B,size,nnz2)
           allocate(Row1(nnz1),Col1(nnz1),Val1(nnz1),Row2(nnz2),Col2(nnz2),Val2(nnz2))
           call dense2sym_sparse(Row1, Col1,Val1, A,size)
           call dense2sym_sparse(Row2, Col2,Val2, B,size)
           if (method.eq.'dense') then 
                      call measure_mult("dense naive", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           elseif(method.eq.'sparse') then 
                      call measure_mult("sparse symmetry", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           else
                      call measure_mult("LAPACK", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           end if
           Nsize(i)= Nmult
           tsize(i)= totalTime
           print*,i 
           deallocate(A,B,C,Row1,Row2,Col1,Col2,Val1,Val2)
           size = size + stepsize 
           i = i + 1 
         end do
        j = i 
        open(1,file = 'sizesparse.dat', status = 'new')
        do i=1,j
          write(1,*) Xsize(i), Nsize(i), tsize(i)
        end do
        close(1)
        i =1 
        stepsize = 1 
        x= int(real(densitymax-densitymin)/stepsize)
        allocate(Xdens(x),tdens(x),Ndens(x))
        size = 100
        density=densitymin
        do while(density < densitymax)
           allocate(A(size,size),B(size,size),C(size,size))
           Xdens(i) = density 
           Nmult = 0 
           totalTime = 0
           call genmat(A,size,density)
           call genmat(B,size,density)
           call get_nnz(A,size,nnz1)
           call get_nnz(B,size,nnz2)
           allocate(Row1(nnz1),Col1(nnz1),Val1(nnz1),Row2(nnz2),Col2(nnz2),Val2(nnz2))
           call dense2sym_sparse(Row1, Col1,Val1, A,size)
           call dense2sym_sparse(Row2, Col2,Val2, B,size)
           if (method.eq.'dense') then
                      call measure_mult("dense naive", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           elseif(method.eq.'sparse') then
                      call measure_mult("sparse symmetry", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           else
                      call measure_mult("LAPACK", A,B,C,size,Nmult,totalTime, Row1 , Row2, Col1, Col2, Val1, Val2, nnz1, &
                      nnz2,debug_flag)
           end if
           Ndens(i)= Nmult
           tdens(i)= totalTime
           deallocate(A,B,C,Row1,Row2,Col1,Col2,Val1,Val2)
           density = density + stepsize
           i = i + 1
         end do
         j = i 
         open(1, file = 'Desitysparse.dat', status ='new')
         do i=1,j 
          write(1, *) Xdens(i),Ndens(i),tdens(i)
         end do
         close(1)
end program scaling 


      



