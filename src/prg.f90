!< author: Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **MSOLV example of use**
!  </span>
!
! @note
!
! The program asks first for the matrix size: ```1``` for a very small matrix, and ```0```
! for a bigger one.The systems are CC HB format systems with elemental matrices.
!
! <h3>Small system -- nnz=18</h3>
!
! The system is (example provided by MUMP):
! $$ A1 = \begin{pmatrix}
!           -1 & 2 & 3 \\
!            2 & 1 & 1 \\
!            1 & 1 & 1
!          \end{pmatrix} ~
!    A2 = \begin{pmatrix}
!            2 & -1 &  3 \\
!            1 &  2 & -1 \\
!            3 &  2 &  1
!         \end{pmatrix}  \\ \text{ } \\
! \rightarrow a\_elt = (-1, 2, 1, 2, 1, 1, 3, 1, 1, 2, 1, 3, -1, 2, 2, 3, -1, 1) \\ \text{ } \\
!  A = \begin{pmatrix}
!           -1 & 2 & 3 & 0 & 0 \\
!            2 & 1 & 1 & 0 & 0 \\
!            1 & 1 & 3 &-1 & 3 \\
!            0 & 0 & 1 & 2 &-1 \\
!            0 & 0 & 3 & 2 & 1
!      \end{pmatrix}  \\ \text{ } \\
! \rightarrow eltvar = (1, 2, 3, 3, 4, 5) \text{ and } eltptr = (1, 4, 7)
! $$
!
! + *eltvar* locates the elemental matrix line in the assembled matrix
! + *eltptr* gives the elemental matrix first entry position in *eltvar*
!          (last position being size(*eltvar*)+1) <br/><br/>
!
! The rhs is (12, 7, 23, 6, 22), and the solution (1, 2, 3, 4, 5)
!
! <h3>Big (well, in fact, medium) system -- nnz=2,097,152</h3>
!
! The system results from a *MUSST* study case and there is of course no theoretical solution to compare with
!
! @endnote
!
! @warning
!
! Some solver implementation are C written, so the arrays may begin at 0. It explains the variable ```dec``` in [[prod_a_x]]
!
! @endwarning
program test_solvers
use omp_lib,   only :   omp_get_wtime
use data_arch, only :   I4, R4, R8, get_unit

use num_param, only :   OPU,          & ! *default output unit*
                        IPU,          & ! *default input unit*
                        SOLV_MESS,    & ! *Solver message: yes=```PRINT_MESS```, no=```NO_MESS```*
                        NO_MESS,      & ! *no output message*
                        PRINT_MESS      ! *solver infos output*

use solver,    only :   MAT_SOLV,     & ! *system data type*
                        solve_syst,   & ! *subroutine for the system resolution*
                        MUMP,         & ! *integer*
                        UMFP,         & ! *integer*
                        SULU,         & ! *integer*
                        convert_matrice_format ! *various matrix format conversion*
implicit none

integer(kind=I4) :: i, ii, uu, size_a_elt, state
real(kind=R8)    :: time1, time2, error
type(MAT_SOLV)   :: slv_struct

real(kind=R8),    dimension(:), allocatable :: b            !! *system right hand side*
real(kind=R8),    dimension(:), allocatable :: a_elt        !! *matrix non-zero entries in CC*
real(kind=R8),    dimension(:), allocatable :: a_elt_ref    !! *initial ```a_elt```*
integer(kind=I4), dimension(:), allocatable :: eltptr       !! *elemental matrices position in ```eltvar```*
integer(kind=I4), dimension(:), allocatable :: eltvar       !! *elemental matrix global lines*

   SOLV_MESS = NO_MESS!PRINT_MESS!
   
   write(OPU, *) 'small matrix? 0:n, 1:y'
   read( IPU, *) i
   
   ! ======================== SOLVER TYPE ===========================================
   do
      write(OPU, *) 'data read, choose solver n° : 0-MA48, 1-SULU, 2-MUMP (ref), 3-UMFP'
      read( IPU, *) ii
      slv_struct%slv_t = ii
      if (ii>=0 .and. ii<=3) exit
   enddo

   ! ======================== SYST DATA READ ===========================================
   call get_unit(uu)
   
   if (i==0) then
      open(uu, file="mat/big_syst.sys")
   else
      open(uu, file="mat/small_syst.sys")
   endif

   ! ======================== SYST INFO ===========================================
   read(uu, *) slv_struct%nn, slv_struct%ne, slv_struct%nvar, slv_struct%nt
   write(OPU, *) '*************** INFO ********************'
   write(OPU, *) 'system size:                  ', slv_struct%nn
   write(OPU, *) 'number of elemental matrices: ', slv_struct%ne
   write(OPU, *) 'number of nnz entries:        ', slv_struct%nt
   write(OPU, *) '*****************************************'
   slv_struct%first = .true.

   allocate( eltvar(slv_struct%nvar  ) )
   allocate( a_elt(     slv_struct%nt    ), stat = state ) ; if (state/=0) stop 'Memory allocation problem in PRG'
   allocate( a_elt_ref( slv_struct%nt    ), stat = state ) ; if (state/=0) stop 'Memory allocation problem in PRG'
   allocate( eltptr(slv_struct%ne +1 ) )
   allocate(      b(slv_struct%nn    ) )

   allocate( slv_struct%eltvar(slv_struct%nvar  ) )
   allocate( slv_struct%a_elt( slv_struct%nt    ), stat = state ) ; if (state/=0) stop 'Memory allocation problem in PRG'
   allocate( slv_struct%eltptr(slv_struct%ne +1 ) )

   call solve_syst(mat = slv_struct, step = 'ini')

   do ii = 1, slv_struct%nvar
      read(uu, *) slv_struct%eltvar(ii)
   enddo
   do ii = 1, slv_struct%ne +1
      read(uu, *) slv_struct%eltptr(ii)
   enddo
   do ii = 1, slv_struct%nt
      read(uu, *) slv_struct%a_elt(ii)
   enddo
   do ii = 1, slv_struct%nn
      read(uu, *) slv_struct%b(ii)
   enddo
   
   close(uu)

   ! The matrices are in CC HB format. Only MUMPS accepts elemental entries, so the
   ! following subroutine converts elemental marices to assembled vectors.
   ! If MUMPS is chosen, nothing is done.
   call convert_matrice_format(mat = slv_struct)
   ! ======================== backup ============================================
   eltvar     = slv_struct%eltvar
   a_elt_ref  = slv_struct%a_elt
   eltptr     = slv_struct%eltptr
   b          = slv_struct%b

   ! ======================== PROCESS ===========================================
   time1 = omp_get_wtime()
      call solve_syst(mat = slv_struct, step = 'ana') ; write(OPU, *) 'system analyzed' ; slv_struct%first = .false.
      call solve_syst(mat = slv_struct, step = 'fac') ; write(OPU, *) 'system factorized'
      call solve_syst(mat = slv_struct, step = 'sol') ; write(OPU, *) 'system solved'
      call solve_syst(mat = slv_struct, step = 'fre') ; write(OPU, *) 'system freed'
   time2 = omp_get_wtime()

   call verif_solution(slv_struct = slv_struct, a_elt = a_elt_ref, b = b, error = error)

   write(OPU,*) 'max error    = ', error
   write(OPU,*) 'elapsed time = ', real(abs(time1 -time2), kind=R4)

   ! Here, the matrix coefficients are modified, but the sparsity is conserved. It gives a look to the ability
   ! of the solver to exploit the symbolic calculations performed before.
   do

      if (slv_struct%slv_t==MUMP) then
         size_a_elt   = size(slv_struct%a_elt)
      else
         size_a_elt   = slv_struct%nt
      endif

      ! the original entries are retrieved, then modified
      a_elt = a_elt_ref
      call modify_a_elt(tab = a_elt, nz = size_a_elt)
      slv_struct%a_elt  = a_elt
      slv_struct%b      = b
      slv_struct%x      = 0._R8

      time1 = omp_get_wtime()
         call solve_syst(mat = slv_struct, step = 'fac')
         call solve_syst(mat = slv_struct, step = 'sol')
         call solve_syst(mat = slv_struct, step = 'fre')
      time2 = omp_get_wtime()

      call verif_solution(slv_struct = slv_struct, a_elt = a_elt, b = b, error = error)

      write(OPU,*) 'max error    = ', error
      write(OPU,*) 'elapsed time = ', real(abs(time1 -time2), kind=R4)

      write(OPU, *) 'STOP? (y=1, n=0)'
      read( IPU, *) ii
      if (ii==1) exit
   enddo

   deallocate(eltvar, a_elt, a_elt_ref, eltptr, b)

   call solve_syst(mat = slv_struct, step = 'end')
   
stop

contains

   !=========================================================================================
   !< @note The product \(\{y\}\) of the system matrix \([A]\) by the solution \(\{x\}\), is
   !        calculated, and compared to the right hand side \(\{b\}\).
   !        The calculated error is the absolute error in %.
   !-----------------------------------------------------------------------------------------
   subroutine verif_solution(slv_struct, a_elt, b, error)
   implicit none
   type(MAT_SOLV), intent(in) :: slv_struct
   real(kind=R8), dimension(:), intent(in) :: a_elt
   real(kind=R8), dimension(:), intent(in) :: b
   real(kind=R8), intent(out) :: error
      real(kind=R8), dimension(slv_struct%nn) :: y
      ! to assess the accuracy, the solution is applied to the
      ! system matrix and compared to the rhs.
      if (slv_struct%slv_t==MUMP) then
         call prod_elemental_x(n      = slv_struct%nn,      &
                               nz     = slv_struct%nt,      &
                               nelt   = slv_struct%ne,      &
                               nvar   = slv_struct%nvar,    &
                               x      = slv_struct%x,       &
                               y      = y,                  &
                               a_elt  = a_elt,              &
                               eltptr = slv_struct%eltptr,  &
                               eltvar = slv_struct%eltvar)
      else
         call prod_a_x(n  = slv_struct%nn,         &
                       nz = slv_struct%nz,         &
                       x  = slv_struct%x,          &
                       y  = y,                     &
                       a_elt = a_elt,              &
                       irow  = slv_struct%irow,    &
                       jptr  = slv_struct%jptr,    &
                       slvt  = slv_struct%slv_t)
      endif

      error = 100*maxval( abs(y(1:slv_struct%nn) -b(1:slv_struct%nn)) )/ &
                  maxval( abs(y(1:slv_struct%nn) +b(1:slv_struct%nn)) )
   return
   endsubroutine verif_solution
   

   !=========================================================================================
   !> @note multiplication of the system coefficient by a random factor
   !-----------------------------------------------------------------------------------------
   subroutine modify_a_elt(tab, nz)
   implicit none
   integer(kind=I4), intent(in) :: nz
   real(kind=R8), intent(inout) :: tab(1:nz)
      real(kind=R8), allocatable :: tmp(:)
      allocate(tmp(1:nz))
      call random_number(harvest=tmp(1:nz))
      tmp(1:nz) = 2*tmp(1:nz) -1.0_R8
      tab(1:nz) = tab(1:nz)*tmp(1:nz)
      deallocate(tmp)
   return
   endsubroutine modify_a_elt

   
   !=========================================================================================
   !> @note \([A] \{x\}\), assembled CC format
   !-----------------------------------------------------------------------------------------
   subroutine prod_a_x(n, nz, x, y, a_elt, irow, jptr, slvt)
   implicit none
   integer(kind=I4), intent(in) :: n, nz, slvt
   real(kind=R8), dimension(nz), intent(in)     :: a_elt
   integer(kind=I4), dimension(nz ), intent(in) :: irow
   integer(kind=I4), dimension(n+1), intent(in) :: jptr
   real(kind=R8), dimension(n), intent(in)      :: x
   real(kind=R8), dimension(n), intent(out)     :: y
      integer(kind=I4) :: i, k, dec
      y(1:n) = 0._R8
      dec = 0
      if (slvt==UMFP .or. slvt==SULU) dec = 1
      do i = 1, n
         do k = jptr(i), jptr(i+1)-1
            y(irow(k +dec) +dec) = y(irow(k +dec) +dec) + x(i)*a_elt(k +dec)
         enddo
      enddo
   return
   endsubroutine prod_a_x

   
   !=========================================================================================
   !> @note \([A] \{x\}\), elemental CC format
   !-----------------------------------------------------------------------------------------
   subroutine prod_elemental_x(n, nz, nelt, nvar, x, y, a_elt, eltptr, eltvar)
   implicit none
   integer(kind=I4), intent(in) :: n, nz, nelt, nvar
   real(kind=R8),    dimension(nz     ), intent(in) :: a_elt
   integer(kind=I4), dimension(nelt +1), intent(in) :: eltptr
   integer(kind=I4), dimension(nvar   ), intent(in) :: eltvar
   real(kind=R8), dimension(n), intent(in)      :: x
   real(kind=R8), dimension(n), intent(out)     :: y
      integer(kind=I4) :: i, j, k, kk, inc_nz, inc_nn, n_elem, i_elem, max_n_elem
      real(kind=R8), dimension(:), allocatable :: a_elt_tmp, x_tmp, y_tmp
      inc_nz = 0
      inc_nn = 0
      n_elem = 0
      y      = 0._R8
      
      max_n_elem = 0
      do i_elem = 1, nelt
         max_n_elem = max(max_n_elem, eltptr(i_elem +1) -eltptr(i_elem))
      enddo

      allocate( a_elt_tmp(1:max_n_elem**2) )
      allocate(     x_tmp(1:max_n_elem   ) )
      allocate(     y_tmp(1:max_n_elem   ) )

      do i_elem = 1, nelt                                               ! browse all elemental matrices
         inc_nn = inc_nn +n_elem                                        ! step in eltvar for matrix number i_elem
         inc_nz = inc_nz +n_elem**2                                     ! step in a_elt for matrix number i_elem
         n_elem = eltptr(i_elem +1) -eltptr(i_elem)                     ! elemental matrix size
         a_elt_tmp(1:n_elem**2) = a_elt(inc_nz +1 : inc_nz +n_elem**2)  ! elemental matrix coefficients
         ! --- elemental rhs
         kk = 0
         do k = inc_nn +1, inc_nn +n_elem
            kk = kk +1
            x_tmp(kk) = x(eltvar(k))
         enddo
         ! --- elemental product
         y_tmp(1:n_elem) = 0._R8
         do i = 1, n_elem
            do j = 1, n_elem
               y_tmp(i) = y_tmp(i) +a_elt_tmp(i +n_elem*(j-1))*x_tmp(j)
            enddo
         enddo
         ! --- elemental product in global vector y
         kk = 0
         do k = inc_nn +1, inc_nn +n_elem
            kk = kk +1
            y(eltvar(k)) = y(eltvar(k)) +y_tmp(kk)
         enddo
      enddo
      deallocate(a_elt_tmp, x_tmp, y_tmp)
   return
   endsubroutine prod_elemental_x

endprogram test_solvers
