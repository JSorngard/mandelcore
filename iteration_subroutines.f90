real*8 function mandel_calc_scaled(cr,ci,maxiters,depth)
!Iterates the mandelbrot function for cr+i*ci until either
!convergence of maxiters iterations.
!Returns a double between 0 and 1.
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

!Pre-check if c is in the main cardioid or the period two bulb.
if((cr +1.d0)**2.d0 + cisqr <= 0.0625d0 .or. mag2*(8.d0*mag2-3.d0) <= .09375d0 - cr) then
    mandel_calc_scaled = 0.d0
    return
end if

!Variable initialization
zr = 0.d0
zi = 0.d0
iters = 0
zrsqr = zr*zr
zisqr = zi*zi

!This loop uses only three multiplications,
!which is the minimum. Continues until |z|^2 <= 6
!instead of 2 because it gives smoother colours.
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

subroutine iterate(re,im,maxiters,depth,result,n,m)
!Results in an array of doubles between 0 and 1 reflecting
!how many iterations of the mandelbrot function are required
!for each point to escape the fractal. 0 if in the set.

!re: an array containing the real parts of the numbers to be iterated.
!im: an array containing the imaginary parts of the numbers to be iterated.
!maxiters: the maximum number of iterations.
!depth: the colour depth that will be used when colouring the image.
!result: on exit contains the result of the iteration.
!n,m: the dimensions of grid. If using f2py you should not need to specify these.

use omp_lib
implicit none
integer,intent(in) :: n,m,maxiters,depth
real*8,intent(in),dimension(n) :: re
real*8,intent(in),dimension(m) :: im
real*8,intent(out),dimension(m,n) :: result

real*8 :: mandel_calc_scaled
integer :: i,j

!If the image is small enough, spreading the work out on multiple cores is unnecesary and makes the computation slower.
if(m < 100) then
    do i=1,n
        do j=1,m
            !(j,i) not (i,j) since python will read the result as transposed.
            result(j,i) = mandel_calc_scaled(re(i),im(j),maxiters,depth)
        end do
    end do
else
    !$OMP parallel do shared(result,re,im)
    do i=1,n
        do j=1,m
            result(j,i) = mandel_calc_scaled(re(i),im(j),maxiters,depth)
        end do
    end do
    !$OMP end parallel do
end if
end subroutine iterate

subroutine iterate_supersampled(re,im,maxiters,depth,samplingfactor,deltar,deltai,result,n,m)
!Same as iterate, but samples every pixel in many
!different places.
!samplingfactor: the number of points that will be samplex along each axis.
!if samplingfactor is 3 a total of 9 points will be sampled.
!deltar: the distance between two pixels along the real axis.
!deltai: the distance between two pixels along the imaginary axis.
use omp_lib
implicit none
integer,intent(in) :: n,m,maxiters,depth,samplingfactor
real*8,intent(in) :: deltar,deltai
real*8,intent(in),dimension(n) :: re
real*8,intent(in),dimension(m) :: im
real*8,intent(out),dimension(m,n) :: result

real*8 :: mandel_calc_scaled,total,coloffset,rowoffset,esc,invfactor
integer :: i,j,k

invfactor = 1.d0/real(samplingfactor,kind=8)

if(m < 100) then
    do i=1,n
        do j=1,m
            total = 0.d0
            do k=1,samplingfactor**2
                !Computes offsets. These should range from -1/samplingfactor
                !to 1/samplingfactor with a 0 included if samplingfator is odd.
                coloffset = (real(mod(k,samplingfactor),kind=8)-1.d0)*invfactor
                rowoffset = (real((k-1)/samplingfactor,kind=8)-1.d0)*invfactor
                esc = mandel_calc_scaled(re(i)+rowoffset*deltar,im(j)+coloffset*deltai,maxiters,depth)
                total = total + esc**2.d0
            end do
            result(j,i) = total/real(samplingfactor**2,kind=8)
        end do
    end do
else    
    !$OMP parallel do shared(result,re,im) private(total,esc,coloffset,rowoffset)
    do i=1,n
        do j=1,m
            total = 0.d0
            do k=1,samplingfactor**2
                coloffset = (real(mod(k,samplingfactor),kind=8)-1.d0)*invfactor
                rowoffset = (real((k-1)/samplingfactor,kind=8)-1.d0)*invfactor
                esc = mandel_calc_scaled(re(i)+rowoffset*deltar,im(j)+coloffset*deltai,maxiters,depth)
                total = total + esc**2.d0
            end do
            result(j,i) = total/real(samplingfactor**2,kind=8)
        end do
    end do
    !$OMP end parallel do
end if


end subroutine iterate_supersampled