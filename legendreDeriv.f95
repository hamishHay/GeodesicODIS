subroutine LegendreDeriv(LegendreP, LegendreDPDX, LMAX, CosColat)
  implicit none

  interface
    subroutine PlmBar_d1(p,dp,lmax,z,csphase,cnorm,exitstatus)
        real*8, intent(out) :: p(:)
        real*8, intent(out) :: dp(:)
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in), optional :: csphase
        integer, intent(in), optional :: cnorm
        integer, intent(out), optional :: exitstatus
    end subroutine PlmBar_d1
  end interface

  integer, intent(in) :: LMAX
  real*8, intent(out) :: LegendreP((LMAX+1)*(LMAX+2)/2)
  real*8, intent(out) :: LegendreDPDX((LMAX+1)*(LMAX+2)/2)
  real*8, intent(in)  :: CosColat
  integer :: returnVAL
  integer :: m,l,l2,m2

  ! write(*,*) LMAX
  ! write(*,*) CosColat


  call PlmBar_d1 (p=LegendreP, dp=LegendreDPDX, lmax=LMAX, z=CosColat, csphase=-1)

  ! write(*,*) 'COS:'

  ! do l=1,LMAX+1
  !  do m=1,l
  !   l2 = l-1
  !   m2 = m-1
  !   write(*,*) l-1,m-1,LegendreDPDX(l2*(l2+1)/2 + m2 + 1)
  !  end do
  ! end do


  return
end subroutine LegendreDeriv
