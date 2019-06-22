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

!$OMP PARALLEL DO shared(grid)
do j=1,m
    do i=1,n
        grid(i,j) = dcmplx(mandel_calc(dble(grid(i,j)),aimag(grid(i,j)),maxiters),0.d0)
    end do
end do
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO shared(grid)
do j=1,m
    do i=1,n
        grid(i,j) = dcmplx(mandel_calc_scaled(dble(grid(i,j)),aimag(grid(i,j)),maxiters,depth),0.d0)
    end do
end do
!$OMP END PARALLEL DO
end subroutine mandel_calc_array_scaled

subroutine fcolour(depth,image,iters,n,m)
use omp_lib
implicit none
integer,intent(in) :: n,m
integer,intent(in) :: depth
real*8,dimension(n,m),intent(in) :: iters
real*8,dimension(n,m,3),intent(inout) :: image
!f2py intent(in,out) :: image
integer :: k,l
real*8 :: T

!$OMP PARALLEL DO shared(image,iters)
do l=1,m
    do k=1,n
        T=iters(k,l)
        !call get_colour(image(k,l,:),iters(k,l),depth)
        image(k,l,1) = T * (depth**(1.d0 - (T**45.d0) * 2.d0))
        image(k,l,2) = (T * 70.d0) - (880.d0 * (T**18.d0)) + (701.d0 * (T**9.d0))
        image(k,l,3) = (T * 80.d0) + ((T**9.d0) * depth) - (950.d0 * (T**99.d0))
    end do
end do
!$OMP END PARALLEL DO


end subroutine fcolour

subroutine get_colour(rgb,scaled_iters,depth)
!Returns an rgb triplet given the result of the
!mandel_calc_scaled function and the colour depth.
implicit none
real*8,intent(in) :: scaled_iters
integer,intent(in) :: depth
real*8,dimension(3),intent(out) :: rgb

rgb(1) = scaled_iters * (depth**(1.d0 - (scaled_iters**45.d0) * 2.d0))
rgb(2) = (scaled_iters * 70.d0) - (880.d0 * (scaled_iters**18.d0)) + (701.d0 * (scaled_iters**9.d0))
rgb(3) = (scaled_iters * 80.d0) + ((scaled_iters**9.d0) * depth) - (950.d0 * (scaled_iters**99.d0))

end subroutine get_colour

!The below routines are based on fasticonv by Sebastian Beyer at https://github.com/sebastianbeyer/fasticonv but modified to work with static memory and openmp.

subroutine naiveGauss (source, filtered, r, nx, ny)
    use omp_lib
    implicit none
    integer, intent(in)                          :: r,nx,ny
    double precision, intent(in)                 :: source(nx,ny)
    double precision, intent(out)                :: filtered(nx,ny)

    double precision, parameter                  :: PI = 4.*atan(1.)

    integer                                      :: i, j, k, l
    integer                                      :: ii, jj      ! inside the array

    double precision                             :: val, wsum
    double precision                             :: dsq         ! distance squared
    double precision                             :: wght        ! weight
    
    !$OMP PARALLEL DO shared(source,filtered)
    do i = 1, nx
        do j = 1, ny
            val = 0
            wsum = 0
            do k = i-r, i+r     ! inner loop over kernel
                do l = j-r, j+r  ! inner loop over kernel
                    ii = min(nx, max(1,k))   ! make sure i is always inside the grid (this implies values are extendet (stretched at the boundaries))
                    jj = min(ny, max(1,l))   ! make sure j is always inside the grid (this implies values are extendet (stretched at the boundaries))
                    dsq = (j-l)**2 + (i-k)**2
                    wght = exp(-dsq / (2.d0*r**2)) / (2.d0*PI*r**2)
                    val = val + source(ii,jj) * wght
                    wsum = wsum + wght
            ! print *, i,j, k, l, ii, jj, dsq
                end do
            end do
            filtered(i,j) = val / wsum
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine naiveGauss


subroutine BoxBlurH (source, filtered, r, nx, ny)
    use omp_lib
    ! computes horizontal blur
    implicit none
    integer, intent(in)                          :: r,nx,ny
    double precision, intent(in)                 :: source(nx,ny)
    double precision, intent(out)                :: filtered(nx,ny)

    double precision                             :: wght  ! weight
    double precision                             :: sum  
    integer                                      :: i,j
    integer                                      :: il    ! leftmost  pixel which should be removed from the accumulator
    integer                                      :: ir    ! rightmost pixel which should be added   to   the accumulator
    
    wght = 1.d0 / (2.d0*r+1.d0)    

    !$OMP PARALLEL DO shared(source,filtered)
    do j = 1, ny   ! loop over all rows
       ! compute sum at first pixel
       sum = source(1,j)
       do i = 1, r  ! loop over box kernel
          sum = sum + source(1,j) + source(i+1,j) ! always take 1 as left pixel, to not get out of grid
          ! print *, sum
       end do

       ! generate output pixel, then update running sum
       do i = 1, nx
          ! print *, j, i, sum
          filtered(i,j) = sum * wght
          il = max(i-r, 1)     ! make sure we dont get off the grid
          ir = min(i+r+1, nx)  ! make sure we dont get off the grid
          sum = sum + source(ir,j) - source(il,j)
       end do
    end do
    !$OMP END PARALLEL DO
       
end subroutine BoxBlurH

  
subroutine BoxBlur (source, filtered, r, nx, ny)
    ! computes box blur, by calling horizontal boxblur twice, the second time with tansposed matrix
    implicit none
    integer, intent(in)                          :: r,nx,ny
    double precision, intent(in)                 :: source(nx,ny)
    double precision, intent(out)                :: filtered(nx,ny)
    
    double precision                             :: source_t(ny,nx)   ! source transposed
    double precision                             :: filtered_t(ny,nx)   ! filtered transposed
 

    ! first horizontal blur
    call BoxBlurH (source, filtered, r, nx ,ny)
    ! then transpose result and call horizontal blur again, now really doing vertical blur
    source_t = transpose (filtered)
    call BoxBlurH (source_t, filtered_t, r, ny, nx)
    ! transpose again, to get back initial shape
    filtered = transpose (filtered_t)
    
end subroutine BoxBlur

  
subroutine fastGauss (source, filtered, r, nx, ny)
    ! computes a fast approximation to a Gaussian filter, by applying a series of box filters
    implicit none
    integer, intent(in)                          :: r,nx,ny
    double precision, intent(in)                 :: source(nx,ny)
    double precision, intent(out)                :: filtered(nx,ny)
    
    double precision                             :: tmpfilter1(nx,ny)   ! tmp array to store intermediate results
    double precision                             :: tmpfilter2(nx,ny)   ! tmp array to store intermediate results
    
    call BoxBlur (source, tmpfilter1, r, nx, ny)
    call BoxBlur (tmpfilter1, tmpfilter2, r, nx, ny) !Possible to use source instead of tmpfilter2 if source has intent(inout), but at the cost of messing up source.
    call BoxBlur (tmpfilter2, filtered, r, nx, ny)

end subroutine fastGauss
