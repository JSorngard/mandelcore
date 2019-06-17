!program test
!implicit none
!real*8 :: mandel_calc
!write(*,*) mandel_calc(-2.d0,0.d0,1000)
!end program test

real*8 function mandel_calc(cr,ci,maxiters)
!An attempt at a fast mandelbrot routine to be used by python through f2py.
implicit none

real*8,intent(in) :: cr,ci
integer,intent(in) :: maxiters
real*8 :: zr,zi,zrsqr,zisqr,q,cisqr
integer :: iters

!Automatically reject points in the main and period-2 bulbs.
cisqr = ci*ci
q = (cr - .25d0)**2.d0 + cisqr
if((cr +1.d0)**2 + cisqr < 0.0625d0 .or. q + cr - .25d0 < .25d0*cisqr) then
    mandel_calc = maxiters
    return
end if

zr = 0
zi = 0
iters = 0
zrsqr = zr*zr
zisqr = zi*zi
do while(zrsqr + zisqr <= 4.d0 .and. iters < maxiters)
    zi = zr*zi
    zi = zi+zi !Multiply by 2.
    zi = zi + ci
    zr = zrsqr - zisqr + cr
    zrsqr = zr*zr
    zisqr = zi*zi
    iters = iters + 1
end do

mandel_calc = iters

end function mandel_calc