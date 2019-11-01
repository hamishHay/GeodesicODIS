subroutine Legendre(LegendreP, LMAX, CosColat)
  implicit none

  interface
    subroutine PlmBar(p,lmax,z,csphase,cnorm,exitstatus)
        real*8, intent(out) :: p(:)
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in), optional :: csphase
        integer, intent(in), optional :: cnorm
        integer, intent(out), optional :: exitstatus
    end subroutine PlmBar
  end interface

  integer, intent(in) :: LMAX
  real*8, intent(out) :: LegendreP((LMAX+1)*(LMAX+2)/2)
  real*8, intent(in)  :: CosColat
  integer :: returnVAL
  integer :: m,l,l2,m2

  call PlmBar (p=LegendreP, lmax=LMAX, z=CosColat, csphase=-1)

  return
end subroutine Legendre
