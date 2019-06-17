!program test
!implicit none
!real*8 :: mandel_calc
!integer,parameter :: n=2
!real*8,dimension(n,n) :: res
!real*8,dimension(n,n,3) :: cres


!write(*,*) mandel_calc(-2.d0,0.d0,1000)

!res(1,1)=0.d0
!res(1,2)=-1.d0
!res(2,1)=-.33d0
!res(2,2)=-.66d0
!cres=0.d0

!call colorize(cres,res,n,n)

!write(*,*) cres(1,1,:)
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

subroutine colorize(colorized,iters,n,m)
implicit none
integer,intent(in) :: n,m
real*8,dimension(n,m),intent(in) :: iters
real*8,dimension(n,m,3),intent(inout) :: colorized
!f2py intent(in,out) :: colorized
integer,parameter :: i=1
integer :: k,l,T

do l=1,m
    do k=1,n
        T=iters(k,l)
        colorized(k,l,1)=T*80.d0 + T**9.d0*i - 950.d0*T**99.d0
        colorized(k,l,2)=T*70.d0 - 880.d0*T**18.d0 + 701.d0*T**9.d0
        colorized(k,l,3)= T*i**(1.d0 - T**45.d0*2.d0)
    end do
end do

end subroutine colorize