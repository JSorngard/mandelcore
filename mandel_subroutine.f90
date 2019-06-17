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

!Automatically reject points in the main cardioid and period-2 bulb.
cisqr = ci*ci
q = (cr - .25d0)**2.d0 + cisqr
if((cr +1.d0)**2.d0 + cisqr < 0.0625d0 .or. q + cr - .25d0 < .25d0*cisqr) then
    mandel_calc = dble(maxiters)
    !mandel_calc = 0.d0
    return
end if

zr = 0.d0
zi = 0.d0
iters = 0
zrsqr = zr*zr
zisqr = zi*zi

!This loop does only 3 multiplies per iteration.
do while(zrsqr + zisqr <= 36.d0 .and. iters < maxiters) !36 instead of 4 for smoother coloring.
    zi = zr*zi
    zi = zi+zi !Multiply by 2.
    zi = zi + ci
    zr = zrsqr - zisqr + cr
    zrsqr = zr*zr
    zisqr = zi*zi
    iters = iters + 1
end do

!if(iters == maxiters) then !In the set
!    mandel_calc = dble(0.d0)
!else
!    mandel_calc = (2+(maxiters-iters)-4*sqrt(zrsqr+zisqr)**(-.4d0))/255.d0 !(0..1]
!end if

mandel_calc = dble(iters)

end function mandel_calc

subroutine colorize(colorized,iters,n,m)
!use omp_lib
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