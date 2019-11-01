subroutine ExtractSHCoeffLL(DATA, N, L_MAX, HARM)
  implicit none

  interface
    subroutine SHExpandDH(grid, n, cilm, lmax, norm, sampling, csphase, lmax_calc, exitstatus)
      real*8, intent(in) :: grid(:,:)
      real*8, intent(out) :: cilm(:,:,:)
      integer, intent(in) :: n
      integer, intent(out) :: lmax
      integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
      integer, intent(out), optional :: exitstatus
    end subroutine SHExpandDH
  end interface

  integer, intent(in) :: N
  real*8, intent(inout) :: DATA(N,2*N)
  integer, intent(in) :: L_MAX
  real*8, intent(out) :: HARM(2,L_MAX+1,L_MAX+1)
  integer :: L_MAX_RETURN
  integer :: returnVAL

  call SHExpandDH (grid=DATA, n=N, cilm=HARM, lmax=L_MAX_RETURN, lmax_calc=L_MAX, csphase=-1, sampling=2, exitstatus=returnVAL)

  return
end subroutine ExtractSHCoeffLL
