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

  call PlmBar_d1 (p=LegendreP, dp=LegendreDPDX, lmax=LMAX, z=CosColat, csphase=-1)

  return
end subroutine LegendreDeriv
