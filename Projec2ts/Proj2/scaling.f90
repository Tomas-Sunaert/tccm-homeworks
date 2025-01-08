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
        logical :: debug_flag
        debug_flag = .false.
        sizemin = 2
        sizemax = 300
        densitymin = 1
        densitymax = 50
        size = sizemin
        density = 25
        method = 'dense'
        i = 1
       
        x= int((sizemax-sizemin)/stepsize)
        allocate(Xsize(x),tsize(x),Nsize(x))
        x= int((densitymax-densitymin)/stepsize)
        allocate(Xdens(x),tdens(x),Ndens(x))

        do while(size < sizemax)
           allocate(A(size,size),B(size,size),C(size,size)) 
           Xdens(i) =  size
          
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
           deallocate(A,B,C,Row1,Row2,Col1,Col2,Val1,Val2)
           size = size + stepsize 
         end do
        size = 100
        density=densitymin
        do while(density < densitymax)
           allocate(A(size,size),B(size,size),C(size,size))
           Xsize(i) =  size

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
         end do

end program scaling 


      



