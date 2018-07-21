!< author: Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **General type definition**
!  </span>
module data_arch
use, intrinsic :: iso_fortran_env, only : output_unit, input_unit, error_unit,   &
                                          int32, int64,              &
                                          real32, real64
implicit none
public ! don't keep global iso_fortran_env parameters

   integer(kind=int32), parameter :: I4 = int32

   integer(kind=I4), parameter :: I8 = int64
   integer(kind=I4), parameter :: R4 = real32
   integer(kind=I4), parameter :: R8 = real64
   integer(kind=I4), parameter :: IN_U  = input_unit
   integer(kind=I4), parameter :: OUT_U = output_unit
   integer(kind=I4), parameter :: ERR_U = error_unit
   integer(kind=I4), parameter :: HIG_I4 = huge(1)

   real(kind=R4), parameter :: PI_R4  = acos(-1._R4)
   real(kind=R8), parameter :: PI_R8  = acos(-1._R8)
   real(kind=R4), parameter :: EPS_R4 = tiny(1._R4)
   real(kind=R8), parameter :: EPS_R8 = tiny(1._R8)
   real(kind=R8), parameter :: HIG_R8 = huge(1._R8)
   real(kind=R8), parameter :: HIG_E8 = log(HIG_R8)
   real(kind=R8), parameter :: EPS_E8 = log(EPS_R8)

   integer(kind=I4), parameter :: EXPO_MAX = exponent(HIG_R8)

   integer(kind=I4) :: NB_THREADS_MAX

private :: output_unit, input_unit, error_unit, int32, int64, real32, real64

contains

!> @note from [John Burkardt website](https://people.sc.fsu.edu/~jburkardt/f_src)
subroutine get_unit(iunit)
implicit none
integer(kind=I4), intent(out) :: iunit
   integer(kind=I4) :: i
   integer(kind=I4) :: ios
   logical(kind=I4) :: lopen
   iunit = 0
   do i = 10, 99

    if (i /= OUT_U .and. i /= IN_U .and. i /= ERR_U) then
      inquire (unit = i, opened = lopen, iostat = ios)
      if (ios == 0) then
         if ( .not. lopen ) then
            iunit = i
            return
         endif
      endif
    endif

   enddo

return
endsubroutine get_unit


endmodule data_arch

