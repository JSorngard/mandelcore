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

cisqr = ci*ci
mag2 = cr*cr+cisqr

!Automatically reject points in the main cardioid and period-2 bulb.
!         reject period 2 bulb                     reject main cardioid.
if((cr +1.d0)**2.d0 + cisqr <= 0.0625d0 .or. mag2*(8.d0*mag2-3.d0) <= .09375d0 - cr) then
    mandel_calc = maxiters
    return
end if
!Alternative rejection:
!q = (cr - .25d0)**2.d0 + cisqr
!if((cr +1.d0)**2.d0 + cisqr <= 0.0625d0 .or. q*(q+(cr-.25d0)) < .25d0*cisqr) then

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

!$OMP parallel do shared(grid)
do j=1,m
    do i=1,n
        grid(i,j) = dcmplx(mandel_calc(dble(grid(i,j)),aimag(grid(i,j)),maxiters),0.d0)
    end do
end do
!$OMP end parallel do

end subroutine mandel_calc_array

real*8 function mandel_calc_scaled(cr,ci,maxiters,depth)
!Same as mandel_calc but returns a double between 0 and 1.
!0 if in the set and increasing to 1 the further out.
!Also scales the computation result based on how far out the point got.
!This results in smoother colour shading.
implicit none
real*8,intent(in) :: cr,ci
integer,intent(in) :: maxiters,depth
real*8 :: zr,zi,zrsqr,zisqr,mag2,cisqr
integer :: iters

cisqr = ci*ci
mag2 = cr*cr+cisqr

if((cr +1.d0)**2.d0 + cisqr <= 0.0625d0 .or. mag2*(8.d0*mag2-3.d0) <= .09375d0 - cr) then
    mandel_calc_scaled = 0.d0 !Using 0 instead of maxiters for points in the set.
    return
end if

zr = 0.d0
zi = 0.d0
iters = 0
zrsqr = zr*zr
zisqr = zi*zi

do while(zrsqr + zisqr <= 36.d0 .and. iters < maxiters)
    zi = zr*zi
    zi = zi + zi
    zi = zi + ci
    zr = zrsqr - zisqr + cr
    zrsqr = zr*zr
    zisqr = zi*zi
    iters = iters + 1
end do

if(iters == maxiters) then !In the set
    mandel_calc_scaled = 0.d0
else
    !This formula scales the answer to a real between 0 and 1 in a way that smoothens out the colours far away.
    mandel_calc_scaled = (2.d0+(maxiters-iters-2.d0)-4.d0*sqrt(zrsqr+zisqr)**(-.4d0))/dble(depth) !Puts the answer in the range 0-1
end if

end function mandel_calc_scaled

subroutine mandel_calc_array_scaled(grid,maxiters,depth,n,m)
!Same as mandel_calc_array but using the scaled version.
!Resultis in an array of doubles instead of integers.
!Stores the result in a grid of complex doubles to save memory space though,
!since complex doubles are needed for the input.
use omp_lib
implicit none
integer,intent(in) :: n,m,maxiters,depth
complex*16,intent(inout),dimension(n,m) :: grid
!f2py intent(in,out) :: grid

real*8 :: mandel_calc_scaled
integer :: i,j

!If the image is small enough, spreading the work out on multiple cores is unnecesary and makes the computation slower.
if(m < 300) then
    do j=1,m
        do i=1,n
            grid(i,j) = dcmplx(mandel_calc_scaled(dble(grid(i,j)),aimag(grid(i,j)),maxiters,depth),0.d0)
        end do
    end do
else
    !$OMP parallel do shared(grid)
    do j=1,m
        do i=1,n
            grid(i,j) = dcmplx(mandel_calc_scaled(dble(grid(i,j)),aimag(grid(i,j)),maxiters,depth),0.d0)
        end do
    end do
    !$OMP end parallel do
end if
end subroutine mandel_calc_array_scaled

subroutine mandel_calc_array_scaled_supersampled(grid,maxiters,depth,samplingfactor,deltar,deltai,n,m)
!Same as mandel_calc_array_scaled, but samples every pixel in many
!different places.
!grid: an array containing the number to be iterated.
!on exit contains the result of the iteration.
!maxiters: the maximum number of iterations.
!depth: the colour depth that will be used when colouring the image.
!samplingfactor: the number of points that will be samplex along each axis.
!if samplingfactor is 3 a total of 9 points will be sampled.
!deltar: the distance between two pixels along the real axis.
!deltai: the distance between two pixels along the imaginary axis.
!n,m: the dimensions of grid. If using f2py you should not need to specify these.
use omp_lib
implicit none
integer,intent(in) :: n,m,maxiters,depth,samplingfactor
real*8,intent(in) :: deltar,deltai
complex*16,intent(inout),dimension(n,m) :: grid
!f2py intent(in,out) :: grid

real*8 :: mandel_calc_scaled,total,coloffset,rowoffset,esc,invfactor
integer :: i,j,k

invfactor = 1.d0/real(samplingfactor,kind=8)

if(m < 300) then
    do j=1,m
        do i=1,n
            total = 0.d0
            do k=1,samplingfactor**2
                !Computes offsets. These should range from -1/samplingfactor
                !to 1/samplingfactor with a 0 included if samplingfator is odd.
                coloffset = (real(mod(k,samplingfactor),kind=8)-1.d0)*invfactor
                rowoffset = (real((k-1)/samplingfactor,kind=8)-1.d0)*invfactor
                esc = mandel_calc_scaled(dble(grid(i,j))+rowoffset*deltar,imag(grid(i,j))+coloffset*deltai,maxiters,depth)
                total = total + esc**2.d0
            end do
            grid(i,j) = total/real(samplingfactor**2,kind=8)
        end do
    end do
else    
    !$OMP parallel do shared(grid) private(total,esc,coloffset,rowoffset)
    do j=1,m
        do i=1,n
            total = 0.d0
            do k=1,samplingfactor**2
                coloffset = (real(mod(k,samplingfactor),kind=8)-1.d0)*invfactor
                rowoffset = (real((k-1)/samplingfactor,kind=8)-1.d0)*invfactor
                esc = mandel_calc_scaled(dble(grid(i,j))+rowoffset*deltar,imag(grid(i,j))+coloffset*deltai,maxiters,depth)
                total = total + esc**2.d0
            end do
            grid(i,j) = total/real(samplingfactor**2,kind=8)
        end do
    end do
    !$OMP end parallel do
end if


end subroutine mandel_calc_array_scaled_supersampled


!impure elemental subroutine pow(base,exponent)
!implicit none
!  !  $    OMP declare simd(pow) uniform(exponent)
!real*8, intent(inout) :: base
!f2py,intent(in,out) :: base
!real*8, intent(in) :: exponent
!base = base**exponent
!end subroutine pow


subroutine multicore_pow(array,exponent,n,m)
!Returns the input array with every element raised to the power of exponent.
!Computes this in multiple threads, so this should only be used on large arrays.
!This subroutine is MUCH slower than numpy.
implicit none
integer,intent(in) :: n,m
real*8,dimension(n,m),intent(inout) :: array
!f2py intent(in,out) :: array
real*8,intent(in) :: exponent
integer :: i,j

!interface
!    impure elemental subroutine pow(base,exponent)
!        !$OMP declare simd(pow) uniform(exponent)
!        real*8,intent(inout) :: base
        !f2py, intent(in,out) :: base
!        real*8,intent(in) :: exponent
!    end subroutine pow
!end interface

!If exponent=1 we don't need to do anything.
if(exponent == 1.d0) then
    return
end if

!$OMP parallel shared(array)
do j=1,m
    do i=1,n
        !call pow(array(n,m),exponent)
        array(n,m) = array(n,m)**exponent
    end do
end do
!$OMP end parallel

end subroutine multicore_pow
