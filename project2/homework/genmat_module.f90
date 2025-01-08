Module genmat_module
        implicit none
        contains

        subroutine genmat(C, size,density)
                integer, intent(in) :: size, density
                real(kind=8),intent(out) :: C(:,:)
                real(kind=8) :: u
                integer :: i,j, num_nonzero, count
                C = 0.0D0
                num_nonzero = int(size*size*density/200.0)
                count =0
                do while (count < num_nonzero)
                    call random_number(u) 
                    i = floor(size*u + 1)
                    call random_number(u)

                    j = floor(size*u + 1)
                    if (C(i,j) == 0.0D0) then
                            call random_number(u)
                            C(i,j) = u
                            C(j,i) = u
                            count = count + 1 
                    end if 
               end do


       end subroutine genmat
end module genmat_module

