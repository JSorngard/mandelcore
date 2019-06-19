!program test
!implicit none
!integer :: mandel_calc,maxiters
!integer,parameter :: n=2
!complex*16,dimension(n,n) :: res

!maxiters=1000

!write(*,*) mandel_calc(.3333d0,0.d0,maxiters)

!res(1,1) = dcmplx(0.d0,0.d0)
!res(1,2) = dcmplx(.333333d0,0.d0)
!res(2,1) = dcmplx(-1.d0,0.d0)
!res(2,2) = dcmplx(0.d0,.5d0)

!call mandel_calc_array(res,maxiters,n,n)

!write(*,*) res

!end program test

integer function mandel_calc(cr,ci,maxiters)
!An attempt at a fast mandelbrot routine to be used by python through f2py.
implicit none

real*8,intent(in) :: cr,ci
integer,intent(in) :: maxiters
real*8 :: zr,zi,zrsqr,zisqr,mag2,cisqr
integer :: iters

!Automatically reject points in the main cardioid and period-2 bulb.
cisqr = ci*ci
mag2 = cr*cr+cisqr

!         reject period 2 bulb                     reject main cardioid.
if((cr +1.d0)**2.d0 + cisqr <= 0.0625d0 .or. mag2*(8.d0*mag2-3.d0) <= .09375d0) then
    mandel_calc = maxiters
    !mandel_calc = 0.d0
    return
end if

!Alternative way of checking main cardioid: compute q
!q = (cr - .25d0)**2.d0 + cisqr
!And check whether
!q*(q+(cr-.25d0)) < .25d0*cisqr
!Might be faster?

zr = 0.d0
zi = 0.d0
iters = 0
zrsqr = zr*zr
zisqr = zi*zi

!This loop does only 3 multiplies per iteration.
do while(zrsqr + zisqr <= 36.d0 .and. iters < maxiters) !36 instead of 4 for smoother coloring.
    zi = zr*zi
    zi = zi + zi !Multiply by 2.
    zi = zi + ci
    zr = zrsqr - zisqr + cr
    zrsqr = zr*zr
    zisqr = zi*zi
    iters = iters + 1
end do

!if(iters == maxiters) then !In the set
!    mandel_calc = dble(0.d0)
!else
!    mandel_calc = (2+(maxiters-iters)-4*sqrt(zrsqr+zisqr)**(-.4d0))/255.d0 !Puts the answer in the range 0-1
!end if

mandel_calc = iters

end function mandel_calc

subroutine mandel_calc_array(grid,maxiters,n,m)
!Multithreaded version of mandel_calc working on the entire array at once.
use omp_lib
implicit none
integer,intent(in) :: n,m,maxiters
complex*16,intent(inout),dimension(n,m) :: grid
!f2py intent(in,out) :: grid

integer :: mandel_calc
integer :: i,j

!$OMP PARALLEL DO shared(grid)
do j=1,m
    do i=1,n
        grid(i,j) = dcmplx(mandel_calc(dble(grid(i,j)),aimag(grid(i,j)),maxiters),0.d0)
    end do
end do
!$OMP END PARALLEL DO

end subroutine mandel_calc_array

subroutine colorize(colorized,iters,n,m)
use omp_lib
implicit none
integer,intent(in) :: n,m
real*8,dimension(n,m),intent(in) :: iters
real*8,dimension(n,m,3),intent(inout) :: colorized
!f2py intent(in,out) :: colorized
integer,parameter :: i=1
integer :: k,l
real*8 :: T

!$OMP PARALLEL DO shared(colorized,iters)
do l=1,m
    do k=1,n
        T=iters(k,l)
        colorized(k,l,1)=T*80.d0 + T**9.d0*i - 950.d0*T**99.d0
        colorized(k,l,2)=T*70.d0 - 880.d0*T**18.d0 + 701.d0*T**9.d0
        colorized(k,l,3)= T*i**(1.d0 - T**45.d0*2.d0)
    end do
end do
!$OMP END PARALLEL DO

end subroutine colorize