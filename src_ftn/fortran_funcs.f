


subroutine laplace(stuff, fuzz, m, n)

implicit none
integer, intent(in) :: m
integer, intent(in) :: n
integer :: j,k
real, intent(in), dimension(m,n) :: stuff
real, intent(out), dimension(m,n) :: fuzz

do j = 1, m -1
  do k = 1, n-1
    fuzz(j,k) = 0.5*stuff(j,k)+stuff(j,k)*stuff(j,k)
  end do
end do
end subroutine laplace


!subroutine bilinear_interp(x_to,y_to,x_from,y_from,var_from,var_to,n,m,l)
!  implicit none
!  integer, intent(in):: n,m,l
!  integer :: i,j,k
!  real, intent(in),dimension(n,m) :: x_to,y_to
!  real, intent(in),dimension(l) :: x_from, y_from, var_from
!  real, intent(out), dimension(n,m) :: var_to
!  do i = 1, m
!    do j=1, n
!      do k=1,l-1
!        fxy1=5.!var_from(k)*(1.-(x_to(i,j)-x_from(k))/(x_from(k+1)-x_from(k)))+var_from(k)*(1.-(x_to(i,j)-x_from(k))/(x_from(k+1)-x_from(k)))
!      enddo
!    enddo
!  enddo
!end subroutine bilinear_interp


subroutine interp1(xData,yData,xVal,yVal,n,m)
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values
  implicit none
  real, intent(in), dimension(n) :: xData, yData
  real, intent(in), dimension(m) :: xVal
  real, intent(out), dimension(m) :: yVal
  integer :: inputIndex, dataIndex
  integer, intent(in) :: n,m
  real :: minXdata, xRange, weight,maxXData

  minXData = xData(1)
  maxXData = xData(n)
  xRange = maxXData - minXData
  do inputIndex = 1, m
      ! possible checks for out of range xVal could go here

      ! this will work if x is uniformly spaced, otherwise increment
      ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
      dataIndex = floor((xVal(inputIndex)-minXData))
      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)

      print *, dataIndex,weight
  end do


end subroutine interp1


subroutine curvature(dem_in,eta,ny,nx,omega_c)
  implicit none
  real, intent(in), dimension(ny,nx) :: dem_in
  real, intent(in) :: eta
  integer, intent(in) :: ny,nx
  real, intent(out), dimension(ny,nx) :: omega_c
  integer :: i,j

  do j=2, ny-1
    do i=2, nx-1
      omega_c(j,i)=1/4.*((dem_in(j,i)-0.5*(dem_in(j-1,i)+dem_in(j+1,i)))/(2.*eta)+  &
        (dem_in(j,i)-0.5*(dem_in(j,i-1)+dem_in(j,i+1)))/(2.*eta)+                   &
        (dem_in(j,i)-0.5*(dem_in(j-1,i-1)+dem_in(j+1,i+1)))/                        &
        (2.*sqrt(2.*eta))+(dem_in(j,i)-0.5*(dem_in(j-1,i+1)+dem_in(j+1,i-1)))/(2.*sqrt(2.*eta)))
    enddo
  enddo

end subroutine curvature

subroutine zinterp(z_in,z_out,v_in,vsfc_in,vele_in,v_out,nx,ny,nz1,nz2)
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values
! This function checks out !

  implicit none
  real, intent(in), dimension(nz1,ny,nx) :: z_in,v_in
  real, dimension(nz1+1,ny,nx) :: z_in1,v_in1
  real, intent(in), dimension(nz2,ny,nx) :: z_out
  real, intent(in), dimension(ny,nx) :: vsfc_in,vele_in
  real, intent(out), dimension(nz2,ny,nx) :: v_out
  integer :: i,j,k,zidx,chk
  integer, intent(in) :: nx,ny,nz1,nz2
  real :: minZ, weight,maxZ,ZZ,slope

  !Perhaps open MP stuff here. !
  do j = 1, ny
    do i = 1, nx
      chk=0
      if (vele_in(j,i) .lt. z_in(1,j,i)) then
        z_in1(1,j,i) = vele_in(j,i)
        v_in1(1,j,i) = vsfc_in(j,i)
        do k = 1, nz1
          z_in1(k+1,j,i) = z_in(k,j,i)
          v_in1(k+1,j,i) = v_in(k,j,i)
        enddo
      elseif (vele_in(j,i) .ge. z_in(nz1,j,i)) then
        do k = 1, nz1
          z_in1(k,j,i) = z_in(k,j,i)
          v_in1(k,j,i) = v_in(k,j,i)
        enddo
        z_in1(nz1+1,j,i) = vele_in(j,i)
        v_in1(nz1+1,j,i) = vsfc_in(j,i)
      else
        do k = 1, nz1-1
          if (vele_in(j,i) .ge. z_in(k,j,i) .and. vele_in(j,i) .lt. z_in(k+1,j,i)) then
            z_in1(k,j,i) = z_in(k,j,i)
            v_in1(k,j,i) = v_in(k,j,i)
            z_in1(k+1,j,i) = vele_in(j,i)
            v_in1(k+1,j,i) = vsfc_in(j,i)
            chk=1
          else
            if (chk .eq. 0) then
              z_in1(k,j,i) = z_in(k,j,i)
              v_in1(k,j,i) = v_in(k,j,i)
            else
              z_in1(k+1,j,i) = z_in(k,j,i)
              v_in1(k+1,j,i) = v_in(k,j,i)
            endif
          endif
        enddo
        z_in1(nz1+1,j,i) = z_in(nz1,j,i)
        v_in1(nz1+1,j,i) = v_in(nz1,j,i)
      endif

  !    print *, z_in1(:,j,i)
  !    print *, z_in(:,j,i)
  !    print *,'-------'
  !    print *, v_in1(:,j,i)
  !    print *, v_in(:,j,i)
  !    print *, vele_in(j,i), vsfc_in(j,i)
  !    stop

      minZ = z_in1(1,j,i)
      maxZ = z_in1(nz1+1,j,i)
      ZZ=minZ
      do k = 1, nz2
      ! possible checks for out of range xVal could go here
      !Set ZZ to Zout at current level index.
        if (z_out(k,j,i) .le. minZ) then !If the z-value is lower than the lowest input Z value at this grid cell, extrapolate down!
          slope=(v_in1(2,j,i)-v_in1(1,j,i))/(z_in1(3,j,i)-z_in1(1,j,i))
          v_out(k,j,i)=v_in1(1,j,i)+slope*(minZ-z_out(k,j,i))

          if(abs(v_out(k,j,i)) .gt. 1000.) then
            print *, "IN lower than min"
            print *,minZ, z_out(k,j,i),slope,v_in1(1,j,i)
          endif

        elseif (z_out(k,j,i) .ge. maxZ) then !If the z-value is higher than the highest input Z value, extrapolate up!
          slope=(v_in1(nz1+1,j,i)-v_in1(nz1,j,i))/(z_in1(nz1+1,j,i)-z_in1(nz1-1,j,i))
          v_out(k,j,i)=v_in1(nz1+1,j,i)+slope*(z_out(k,j,i)-maxZ)

          if(abs(v_out(k,j,i)) .gt. 1000.) then
            print *, "IN Greater than max"
            print *,maxZ, z_out(k,j,i),slope,v_in1(1,j,i)
          endif

        else !It's in range!
          do zidx = 1, nz1
            if (z_out(k,j,i) .ge. z_in1(zidx,j,i) .and. z_out(k,j,i) .lt. z_in1(zidx+1,j,i)) exit
          end do
          ! Found z index
          weight = (z_out(k,j,i) - z_in1(zidx,j,i))/(z_in1(zidx+1,j,i)-z_in1(zidx,j,i))

          v_out(k,j,i) = (1.0-weight)*v_in1(zidx,j,i) + weight*v_in1(zidx+1,j,i)

        endif
      enddo
    enddo
  enddo

end subroutine zinterp
