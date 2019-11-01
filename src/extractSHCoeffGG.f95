subroutine ExtractSHCoeffGG(DATA, LATS, LONS, N_MAX, L_MAX, HARM)
  implicit none

  interface
    subroutine SHExpandLSQ(cilm, d, lat, lon, nmax, lmax, norm, chi2, csphase, exitstatus)
      real*8, intent(in) :: d(:)
      real*8, intent(out) :: cilm(:,:,:)
      real*8, intent(in) :: lat(:), lon(:)
      integer, intent(in) :: nmax
      integer, intent(in) :: lmax
      integer, intent(in), optional :: norm, csphase, chi2
      integer, intent(out), optional :: exitstatus
    end subroutine SHExpandLSQ
  end interface

  integer, intent(in) :: N_MAX
  real*8, intent(in) :: DATA(N_MAX)
  real*8, intent(in) :: LATS(N_MAX)
  real*8, intent(in) :: LONS(N_MAX)
  integer, intent(in) :: L_MAX
  real*8, intent(out) :: HARM(2,L_MAX+1,L_MAX+1)
  integer :: l,m
  integer :: returnVAL

  ! do l=1, N_MAX+1
  !   print *, DATA(l)
  ! enddo

  call SHExpandLSQ (cilm=HARM, d=DATA, lat=LATS, lon=LONS, nmax=N_MAX, lmax=L_MAX, csphase=-1, exitstatus=returnVAL)

  ! do l=1, L_MAX+1
  !   do m=1, l
  !     print *, l-1, m-1, HARM(1, l, m), HARM(2,l,m)
  !   enddo
  ! enddo

  return
end subroutine ExtractSHCoeffGG
