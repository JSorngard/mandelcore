subroutine fcolour(depth,image,iters,n,m)
!Returns an array of RGB triplets given an array of results
!from the scaled algorithms in mandel_subroutine.f90.
!depth: the colour depth of the image.
!This should have been used during the computation.
!iters: an array of containing the result from the mandelbrot iterations.
!image: on exit will contain an array of RGB triplets.
use omp_lib
implicit none
integer,intent(in) :: n,m
integer,intent(in) :: depth
real*8,dimension(n,m),intent(in) :: iters
real*8,dimension(n,m,3),intent(inout) :: image
!f2py intent(in,out) :: image
integer :: k,l
real*8 :: T

!$OMP PARALLEL DO shared(image,iters) private(T)
do l=1,m
    do k=1,n
        T=iters(k,l)

        !Maybe this is faster?
        if(T == 0.d0) then
            image(k,l,:) = 0.d0
        end if
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
!I can not get it to cooperate well with OpenMP with the knowledge I have.
implicit none
real*8,intent(in) :: scaled_iters
integer,intent(in) :: depth
real*8,dimension(3),intent(out) :: rgb

rgb(1) = scaled_iters * (depth**(1.d0 - (scaled_iters**45.d0) * 2.d0))
rgb(2) = (scaled_iters * 70.d0) - (880.d0 * (scaled_iters**18.d0)) + (701.d0 * (scaled_iters**9.d0))
rgb(3) = (scaled_iters * 80.d0) + ((scaled_iters**9.d0) * depth) - (950.d0 * (scaled_iters**99.d0))

end subroutine get_colour

!The below routines are based on fasticonv by Sebastian Beyer at https://github.com/sebastianbeyer/fasticonv but modified to work with static memory and openmp.
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
