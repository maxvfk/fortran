program sub2
use omp_lib    
implicit none
    
    integer n
!$ print *, 'hahahahahahaha'
    n=2
    print *, "z = ",n
read *
end program sub2

subroutine f(x)
    implicit none
    real(kind=8) :: x
    x = x**2
end subroutine f

subroutine g(x)
    implicit none
    real(kind=8) :: x
    x = x+x
end subroutine g
