!< author: Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **Api for different sparse matrix solvers**
!  </span>

module solver
use iso_c_binding, only : C_PTR, C_NULL_PTR, C_ASSOCIATED
use data_arch,     only : I4, I8, R8, EPS_R8, HIG_R8
use num_param,     only : OPU, SOLV_MESS, NO_MESS, PRINT_MESS
use sort_arrays,   only : sort_int_1real, sort_int_1int_1real
!-------------------------------------------------
use hsl_ma48_double
!-------------------------------------------------
use sulu_wrapper,  MAT_SULU => SULU_ENV      , &
           SULU_SAMEPATTERN => SAMEPATTERN   , &
              SULU_FACTORED => FACTORED      , &
                SULU_DOFACT => DOFACT
!-------------------------------------------------
use mumps_wrapper, MAT_MUMP => DMUMPS_STRUC
!-------------------------------------------------
use mumfpack, only : UMFPACK_CONTROL               , &
                     UMFPACK_STATUS                , &
                     UMFPACK_INFO                  , &
                     UMFPACK_PRL                   , &
                     UMFPACK_IRSTEP                , &
                     UMFPACK_NUMERIC_SIZE_ESTIMATE , &
                     UMFPACK_PEAK_MEMORY_ESTIMATE  , &
                     UMFPACK_SIZE_OF_UNIT          , &
                     UMFPACK_DI_DEFAULTS           , &
                     UMFPACK_DI_REPORT_CONTROL     , &
                     UMFPACK_DI_FREE_NUMERIC       , &
                     UMFPACK_DI_FREE_SYMBOLIC      , &
                     UMFPACK_DI_REPORT_INFO        , &
                     S_UMFPACK_DI_SYMBOLIC         , &
                     S_UMFPACK_DI_NUMERIC          , &
                     S_UMFPACK_DI_SOLVE            , &
                     UMFPACK_A
!-------------------------------------------------
implicit none

type MAT_MA48
!! <span style="color:green">All the stuff needed by *HSL_MA48*</span>
   type(ZD11_TYPE)    :: zmat
   type(MA48_CONTROL) :: ctrl
   type(MA48_AINFO)   :: ainf
   type(MA48_FINFO)   :: finf
   type(MA48_SINFO)   :: sinf
   type(MA48_FACTORS) :: fact
   integer(kind=I4)   :: fast
   real(kind=R8), dimension(2):: resid
endtype MAT_MA48

type MAT_UMFP
!! <span style="color:green">All the stuff needed by *UMFPACK*</span>
   type(c_ptr) :: c_symbolic
   type(c_ptr) :: c_numeric
   real(kind=R8), dimension(0:UMFPACK_CONTROL-1) :: c_control
   real(kind=R8), dimension(0:UMFPACK_INFO   -1) :: c_info
endtype MAT_UMFP

type MAT_SOLV
!! <span style="color:green">MUSST high level system type</span>
   integer(kind=I4) :: slv_t                             !! *solver type*
   logical(kind=I4) :: first = .true.                    !! *analysis of the system to be done?*
   integer(kind=I4) :: nn                                !! *number of nodes*
   integer(kind=I4) :: ne                                !! *number of elements*
   integer(kind=I4) :: nt                                !! *number of **a priori** non-zero terms in the matrix*
   integer(kind=I4) :: nz                                !! *number of non-zero terms in the matrix*
   integer(kind=I4) :: nvar                              !! *eltvar length ( if 4-nodes elt -> 2 lines X number of elemental matrices)*
   integer(kind=I4) :: code                              !! *error code*   [not used yet]
   real(kind=R8)    :: error                             !! *error value*  [not used yet]
   character(len=1024) :: mess                           !! *message*      [not used yet]
   !.....................................................
   type(MAT_MUMP) :: matmump                             !! *matrices for mumps solver*
   type(MAT_MA48) :: matma48                             !! *matrices for ma48 solver*
   type(MAT_SULU) :: matsulu                             !! *matrices for SuperLu solver*
   type(MAT_UMFP) :: matumfp                             !! *matrices for Umfpack solver*
   !.....................................................
   integer(kind=I4), dimension(:), allocatable :: eltvar !! *rows in assembled matrix*
   integer(kind=I4), dimension(:), allocatable :: eltptr !! *element rows pointer*
   real(kind=R8),    dimension(:), allocatable :: a_elt  !! *unassembled rigidity matrix*
   !.....................................................
   integer(kind=I4), dimension(:), allocatable :: irow   !! *line number*
   integer(kind=I4), dimension(:), allocatable :: jcol   !! *column number*
   integer(kind=I4), dimension(:), allocatable :: jptr   !! *line pointer*
   !.....................................................
   real(kind=R8), dimension(:), allocatable    :: b      !! *right hand side vector*
   real(kind=R8), dimension(:), allocatable    :: x      !! *unknwon vector*
endtype MAT_SOLV


!< <span style="color:green">MUSST multiscale high level solver type</span>
!  @note ```MS_MAT_SOLV``` is needed by MUSST, but it is useless for the present module
!
type MS_MAT_SOLV
   type(MAT_SOLV) :: ts_mat                                       !! *top-scale solver type matrices*
   type(MAT_SOLV),   dimension(:), allocatable :: bs_mat          !! *bottom-scale solver type matrices (table)*
   integer(kind=I4), dimension(:), allocatable :: ass_loc_in_mat  !! *table for assembly location (for parallel computation)*
endtype ms_mat_solv


! Solver types
integer(kind=I4), parameter :: MA48  = 0 !! *code for     Ma48 solver type*
integer(kind=I4), parameter :: SULU  = 1 !! *code for SUPER LU solver type*
integer(kind=I4), parameter :: MUMP  = 2 !! *code for    MUMPS solver type*
integer(kind=I4), parameter :: UMFP  = 3 !! *code for  UMFPACK solver type*

! What solver to use for bottom- or top- scale grids
integer(kind=I4) :: SOLVER_BS = -1 !! *solver used for bottom scale grids* [not used by the present module]
integer(kind=I4) :: SOLVER_TS = -1 !! *solver used for top scale grids* [not used by the present module]

contains

   !=========================================================================================
   !< @note General hat subroutine that handles the resolution steps:      <br/>
   !  * ```ini``` solver initialization                                    <br/>
   !  * ```ana``` solver analyzis when it's proposed by the solver         <br/>
   !  * ```fac``` solver factorization                                     <br/>
   !  * ```sol``` solver solution                                          <br/>
   !  * ```fre``` solver memory release when it's proposed by the solver   <br/>
   !  * ```end``` solver end
   !-----------------------------------------------------------------------------------------
   subroutine solve_syst(mat, step)
   implicit none
   type(MAT_SOLV), intent(inout) :: mat   !! *high level system type*
   character(len=*), intent(in)  :: step  !! *'ini'=initialize, 'ana'=analyze, 'fac'=factorize, 'sol'=solve, 'fre'=free memory, 'end'=close solver*
      if ( index(step, 'ini')/=0 ) then ; call init_solver(     mat) ; return ; endif
      if ( index(step, 'ana')/=0 ) then ; call analyse_solver(  mat) ; return ; endif
      if ( index(step, 'fac')/=0 ) then ; call factorize_solver(mat) ; return ; endif
      if ( index(step, 'sol')/=0 ) then ; call solution_solver( mat) ; return ; endif
      if ( index(step, 'fre')/=0 ) then ; call freefact_solver( mat) ; return ; endif
      if ( index(step, 'end')/=0 ) then ; call close_solver( mat)    ; return ; endif
      stop 'Bad step chosen in SOLVE_SYST'
   return
   endsubroutine solve_syst


   !=========================================================================================
   !> @note Subroutine to initialize the matrices of the solver
   !-----------------------------------------------------------------------------------------
   subroutine init_solver(mat)
   implicit none
   type(MAT_SOLV), intent(inout) :: mat   !! *high level system type*
      integer(kind=I4) :: ierr

      ! allocation of the system vectors: rhs and unknown
      allocate( mat%b(mat%nn), mat%x(mat%nn) )
      mat%b = HIG_R8
      mat%x = HIG_R8

      ! check solver type
      select case (mat%slv_t)
         case(MA48)
            call ma48_initialize(factors = mat%matma48%fact, &
                                 control = mat%matma48%ctrl)
            select case(SOLV_MESS)
               case(NO_MESS)
                  mat%matma48%ctrl%wp    = -1
                  mat%matma48%ctrl%mp    = -1
                  mat%matma48%ctrl%ldiag = -1
               case(PRINT_MESS:)
                  mat%matma48%ctrl%ldiag = +2
            endselect

         case(MUMP)
            call mpi_init(ierr)
            mat%matmump%comm = MPI_COMM_WORLD
            mat%matmump%job  = -1               ! initialisation
            mat%matmump%sym  =  0               ! no symetry
            mat%matmump%par  =  1               ! MPI, host working
            call dmumps(mat%matmump)
            mat%matmump%icntl(5) =  1           ! Specify element entry : elemental matrices
            if (mat%matmump%infog(1) < 0) then
               write(OPU,'(a,a,i6,a,i9)') ' error return: ',                              &
                                          '  mumps_par%infog(1)= ', mat%matmump%infog(1), &
                                          '  mumps_par%infog(2)= ', mat%matmump%infog(2)
               call mpi_finalize(ierr)
               stop 'MUMP error, INIT_SOLVER'
            endif
            select case(SOLV_MESS)
               case(NO_MESS)
                  mat%matmump%icntl(4) =  1 ! error output only
               case(PRINT_MESS:)
                  mat%matmump%icntl(4) =  3 ! all error output
            endselect

         case(SULU)
            call init_superlu(sulu = mat%matsulu)
            select case(SOLV_MESS)
               case(NO_MESS)
                  mat%matsulu%options%PrintStat = 0
               case(PRINT_MESS:)
                  mat%matsulu%options%PrintStat = 1
            endselect

         case(UMFP)
            mat%matumfp%c_numeric = C_NULL_PTR
            call umfpack_di_defaults(mat%matumfp%c_control)
            select case(SOLV_MESS)
               case(NO_MESS)
                  mat%matumfp%c_control(UMFPACK_PRL) = 1
               case(PRINT_MESS:)
                  mat%matumfp%c_control(UMFPACK_PRL) = 2
            endselect
            call umfpack_di_report_control(mat%matumfp%c_control)

         case default
            stop 'Unknown solver type, INIT_SOLVER'

      endselect

   return
   endsubroutine init_solver


   !=========================================================================================
   !> @note Subroutine to analyse, factorize (symbolic) the matrix of the system
   !-----------------------------------------------------------------------------------------
   subroutine analyse_solver(mat)
   implicit none
   type(MAT_SOLV), intent(inout), target :: mat   !! *high level system type*

      select case(mat%slv_t)

         case(MA48)
            mat%matma48%zmat%row => mat%irow
            mat%matma48%zmat%col => mat%jcol
            mat%matma48%zmat%val => mat%a_elt

            mat%matma48%zmat%n   = mat%nn
            mat%matma48%zmat%m   = mat%nn
            mat%matma48%zmat%ne  = mat%nz

            call ma48_analyse(matrix = mat%matma48%zmat, &
                             factors = mat%matma48%fact, &
                             control = mat%matma48%ctrl, &
                               ainfo = mat%matma48%ainf, &
                               finfo = mat%matma48%finf)
            if (mat%matma48%ainf%flag < 0) then
               write(OPU,*) 'Failure of ma48_analyse with ainfop%flag = ', mat%matma48%ainf%flag
               stop
            endif

         case(MUMP)
            mat%matmump%eltptr   => mat%eltptr
            mat%matmump%eltvar   => mat%eltvar
            mat%matmump%a_elt    => mat%a_elt
            mat%matmump%rhs      => mat%b

            mat%matmump%n    = mat%nn
            mat%matmump%nelt = mat%ne

            mat%matmump%job = 1 ! performs the analysis

            call dmumps(mat%matmump)

         case(SULU)
            mat%matsulu%irow   => mat%irow
            mat%matsulu%jptr   => mat%jptr
            mat%matsulu%a_elt  => mat%a_elt
            mat%matsulu%b      => mat%b

            mat%matsulu%n     = mat%nn
            mat%matsulu%nz    = mat%nz
            mat%matsulu%nrhs  = 1
            mat%matsulu%first = .true.

            call prep_superlu(sulu = mat%matsulu)

         case(UMFP)
            call s_umfpack_di_symbolic(n_row    = mat%nn,                  &
                                       n_col    = mat%nn,                  &
                                       Ap       = mat%jptr,                &
                                       Ai       = mat%irow,                &
                                       Ax       = mat%a_elt,               &
                                       Symbolic = mat%matumfp%c_symbolic,  &
                                       Control  = mat%matumfp%c_control,   &
                                       Info     = mat%matumfp%c_info)

         case default
            stop 'ANALYSE_SOLVER : unknown solver type'

      endselect

   return
   endsubroutine analyse_solver


   !=========================================================================================
   !> @note Subroutine to factorize the matrix of the system
   !-----------------------------------------------------------------------------------------
   subroutine factorize_solver(mat)
   implicit none
   type(MAT_SOLV), intent(inout) :: mat   !! *high level system type*

      select case(mat%slv_t)

         case(MA48)
            call ma48_factorize( matrix  = mat%matma48%zmat, &
                                 factors = mat%matma48%fact, &
                                 control = mat%matma48%ctrl, &
                                 finfo   = mat%matma48%finf, &
                                 fast    = mat%matma48%fast)
            if (mat%matma48%finf%flag < 0) then
               write(OPU,*) 'Failure of ma48_factorize with finfo%flag = ', mat%matma48%finf%flag
               stop
            endif

         case(MUMP)
            mat%matmump%job = 2
            call dmumps(mat%matmump)

         case(SULU)
            ! Just factorize once, in the case that the matrix doesn't change much
            if (mat%matsulu%first) call fact_superlu( sulu    = mat%matsulu, &
                                                      verbose = (mat%matsulu%options%PrintStat==1))

         case(UMFP)
            call umfpack_di_free_numeric(Numeric = mat%matumfp%c_numeric) ! first release memory
            call s_umfpack_di_numeric(Ap       = mat%jptr,                 &
                                      Ai       = mat%irow,                 &
                                      Ax       = mat%a_elt,                &
                                      Symbolic = mat%matumfp%c_symbolic,   &
                                      Numeric  = mat%matumfp%c_numeric,    &
                                      Control  = mat%matumfp%c_control,    &
                                      Info     = mat%matumfp%c_info)

         case default
            stop 'Unknown solver type, FACTORIZE_SOLVER'

      endselect
   return
   endsubroutine factorize_solver


   !=========================================================================================
   !> @note Subroutine to solve the system \([A]\{x\} = \{b\}\) (sparse A)
   !-----------------------------------------------------------------------------------------
   subroutine solution_solver(mat)
   implicit none
   type(MAT_SOLV), target, intent(inout) :: mat   !! *high level system type*
      integer(kind=I4) :: ierr

      select case(mat%slv_t)

         case(MA48)
            call ma48_solve(matrix = mat%matma48%zmat,   &
                           factors = mat%matma48%fact,   &
                               rhs = mat%b,              &
                                 x = mat%x,              &
                           control = mat%matma48%ctrl,   &
                             sinfo = mat%matma48%sinf,   &
                             resid = mat%matma48%resid,  &
                             error = mat%error)

         case(MUMP)
            mat%matmump%job = 3
            call dmumps(mat%matmump)
            if (mat%matmump%infog(1) < 0) then
               write(OPU,'(a,a,i6,a,i9)') ' error return: ',                              &
                                          '  mumps_par%infog(1)= ', mat%matmump%infog(1), &
                                          '  mumps_par%infog(2)= ', mat%matmump%infog(2)
               call mpi_finalize(ierr)
               stop 'Error in SOLUTION_SOLVER'
            endif
            ! solution has been assembled on the host
            if ( mat%matmump%myid == 0 ) then
               mat%x(1:mat%nn) = mat%matmump%rhs(1:mat%nn)
            endif

         case(SULU)
            call solv_superlu(sol_x = mat%x,        &
                              sulu  = mat%matsulu,  &
                           verbose  = (mat%matsulu%options%PrintStat==1))
            mat%matsulu%first = .false.

         case(UMFP)
            ! Numeric factors must exist
            if (.not.C_ASSOCIATED(mat%matumfp%c_numeric)) call solve_syst(mat, 'fac')
            mat%matumfp%c_control(UMFPACK_IRSTEP) = 0 ! solve ax=b, without iterative refinement

            ! If you want to evaluate the required RAM (Go)
            ! write(*,*) mat%matumfp%c_info(UMFPACK_PEAK_MEMORY_ESTIMATE)/mat%matumfp%c_info(UMFPACK_SIZE_OF_UNIT)/1e9
            ! write(*,*) sizeof(mat%a_elt)/1e9

            call s_umfpack_di_solve(sys = UMFPACK_A,              &
                                      x = mat%x,                  &
                                      b = mat%b,                  &
                                numeric = mat%matumfp%c_numeric,  &
                                control = mat%matumfp%c_control,  &
                                   info = mat%matumfp%c_info)

            if (mat%matumfp%c_info(UMFPACK_STATUS) < 0) then
               write(OPU, *) 'error occurred in umfpack_di_solve: ', mat%matumfp%c_info(UMFPACK_STATUS)
               stop 'Error in SOLUTION_SOLVER'
            endif

         case default
            stop 'Unknown solver type, SOLUTION_SOLVER'

      endselect

   return
   endsubroutine solution_solver


   !=========================================================================================
   !> @note Subroutine to free the factors if applicable
   !-----------------------------------------------------------------------------------------
   subroutine freefact_solver(mat)
   implicit none
   type(MAT_SOLV), intent(inout) :: mat   !! *high level system type*

      select case(mat%slv_t)

         case(MA48)
            continue

         case(MUMP)
            continue

         case(SULU)
            call free_superlu()

         case(UMFP)
            call umfpack_di_free_numeric(Numeric = mat%matumfp%c_numeric) ! free the numeric factorization

         case default
            stop 'unknown solver type, FREEFACT_SOLVER'

      endselect

   return
   endsubroutine freefact_solver


   !=========================================================================================
   !> @note Subroutine to close the solver
   !-----------------------------------------------------------------------------------------
   subroutine close_solver(mat)
   implicit none
   type(MAT_SOLV), intent(inout) :: mat   !! *high level system type*
      integer(kind=I4) :: ierr

      deallocate( mat%eltptr )
      deallocate( mat%eltvar )

      select case(mat%slv_t)

         case(MA48)
            nullify( mat%matma48%zmat%row )
            nullify( mat%matma48%zmat%col )
            nullify( mat%matma48%zmat%val )

            call ma48_finalize(factors = mat%matma48%fact, &
                               control = mat%matma48%ctrl, &
                               info    = ierr)
            deallocate( mat%b, mat%x )

         case(MUMP)
            if ( mat%matmump%myid == 0 )then
               nullify( mat%matmump%rhs    )
               nullify( mat%matmump%eltptr )
               nullify( mat%matmump%eltvar )
               nullify( mat%matmump%a_elt  )
            endif

            ! destroy the instance (deallocate internal data structures)
            mat%matmump%job = -2
            call dmumps(mat%matmump)
            if (mat%matmump%infog(1) < 0) then
               write(OPU,'(a,a,i6,a,i9)') ' error return: ',                              &
                                          '  mumps_par%infog(1)= ', mat%matmump%infog(1), &
                                          '  mumps_par%infog(2)= ', mat%matmump%infog(2)
               call mpi_finalize(ierr)
               stop 'Error MUMP in close_solver'
            endif
            call mpi_finalize(ierr)
            deallocate( mat%b, mat%x )

         case(SULU)
            call close_superlu(sulu = mat%matsulu)

            nullify( mat%matsulu%b     )
            nullify( mat%matsulu%irow  )
            nullify( mat%matsulu%jptr  )
            nullify( mat%matsulu%a_elt )

         case(UMFP)
            deallocate( mat%jptr, mat%irow )

            call umfpack_di_free_numeric(Numeric = mat%matumfp%c_numeric)     ! free the numeric factorization
            ! no lu factors (numeric) are in memory at this point.
            call umfpack_di_free_symbolic(Symbolic = mat%matumfp%c_symbolic)  ! free the symbolic analysis
            call umfpack_di_report_info(Control = mat%matumfp%c_control,      &
                                        Info    = mat%matumfp%c_info)         ! print final statistics
            deallocate( mat%b, mat%x )

         case default
            stop 'Unknown solver type, close_solver'

      endselect

   return
   endsubroutine close_solver


   !=========================================================================================
   !> @note Subroutine to transform the Rutherford Boeing format into Harwell Boeing and triplet
   !-----------------------------------------------------------------------------------------
   subroutine convert_matrice_format(mat)
   implicit none
   type(MAT_SOLV), intent(inout), target :: mat   !! *high level system type*
      integer(kind=I4) :: i

      ! =======================================================================================================================
      ! Compressed Column Storage (CCS) is also called the Harwell-Boeing sparse matrix format.
      ! ************************************************
      ! Elemental entries (example provided by MUMP):
      ! A1 = 1|-1  2  3| A2 = 3|2 -1  3|
      !      2| 2  1  1|      4|1  2 -1|
      !      3| 1  1  1|      5|3  2  1| => a_elt = (-1, 2, 1, 2, 1, 1, 3, 1, 1, 2, 1, 3, -1, 2, 2, 3, -1, 1)
      !
      ! A  = 1|-1  2  3  0  0|
      !      2| 2  1  1  0  0|
      !      3| 1  1  3 -1  3|
      !      4| 0  0  1  2 -1|
      !      5| 0  0  3  2  1| => eltvar = (1, 2, 3, 3, 4, 5), it locates the elemental matrix line in the assembled matrix
      !                        => eltptr = (1, 4, 7), it gives the elemental matrix first entry position in eltvar (last
      !                                               position being size(eltvar)+1)
      !
      ! ************************************************
      ! Assembled matrix :
      ! A being the same, a_elt = (-1, 2, 1, 2, 1, 1, 3, 1, 3, 1, 3, -1, 2, 2, 3, -1, 1)
      !                    irow = ( 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 5,  3, 4, 5, 3,  4, 5)
      !                    jptr = (1, 4, 7, 12, 15, 18)
      !
      ! =======================================================================================================================
      ! Triplet form
      ! ************************************************
      ! For each non zero a_elt entry, returns its row and column number
      ! A being the same, a_elt = (-1, 2, 1, 2, 1, 1, 3, 1, 3, 1, 3, -1, 2, 2, 3, -1, 1)
      !                    irow = ( 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 5,  3, 4, 5, 3,  4, 5)
      !                    jcol = ( 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3,  4, 4, 4, 5,  5, 5)

      call from_elemental_to_assembled(mat = mat)

      select case(mat%slv_t)

         case(MA48)  ! Triplet form: irow, jcol, a_elt
            continue

         case(MUMP)  ! Compressed row, elemental entries chosen: eltptr, eltvar, a_elt
            continue

         case(UMFP, SULU)  ! Compressed row, assembled form: irow, jptr, a_elt
            deallocate( mat%jcol ) ! no more needed
            ! UMFP and SULU allocations begin at 0 (C convention)
            do i = 1, mat%nn +1
               mat%jptr(i) = mat%jptr(i) -1
            enddo

            if (mat%slv_t==SULU) mat%matsulu%jptr => mat%jptr ! otherwise, matsulu%jptr will be associated to
                                                              ! a deallocated part of memory

            do i = 1, mat%nz
               mat%irow(i) = mat%irow(i) -1
            enddo

         case default
            stop 'Unknown solver type, CONVERT_MATRICE_FORMAT'

      endselect
   return
   endsubroutine convert_matrice_format


   !=========================================================================================
   !> @note Subroutine to transform the elemental entries into assembled CC vectors
   !-----------------------------------------------------------------------------------------
   subroutine from_elemental_to_assembled(mat)
   implicit none
   type(MAT_SOLV), intent(inout), target :: mat   !! *high level system type*
      integer(kind=I4), pointer :: solver, nb_elt, n, ntot, nz
      integer(kind=I4), dimension(:), pointer :: eltptr
      integer(kind=I4), dimension(:), pointer :: eltvar
      real(kind=R8),    dimension(:), pointer :: a_elt
      integer(kind=I4), dimension(:), pointer :: irow, jcol
      integer(kind=I4), dimension(:), pointer :: jptr
      integer(kind=I4) :: inelt, imatorder, i, j, ii, jj, i1, i2, ir1, ir2, jr1, jr2, itmp, innz, state

      solver   => mat%slv_t
      nb_elt   => mat%ne
      n        => mat%nn
      ntot     => mat%nt
      nz       => mat%nz

      if (solver==MUMP) return

      ! conversion from elemental form to triplet form, perhaps with null a_elts
      !--------------------------------------------------------------------------
      state = 0
      if (.not.allocated(mat%irow)) allocate( mat%irow(ntot), stat = state )
      if (state/=0) stop 'Memory allocation problem, FROM_ELEMENTAL_TO_ASSEMBLED'
      if (.not.allocated(mat%jcol)) allocate( mat%jcol(ntot), stat = state )
      if (state/=0) stop 'Memory allocation problem, FROM_ELEMENTAL_TO_ASSEMBLED'

      eltptr   => mat%eltptr
      eltvar   => mat%eltvar
      a_elt    => mat%a_elt
      irow     => mat%irow
      jcol     => mat%jcol
      jptr     => mat%jptr

      irow(1:ntot) = -1
      jcol(1:ntot) = -1

      ii = 1
      do inelt = 1, nb_elt
         imatorder = eltptr(inelt +1) -eltptr(inelt)
         do i = 1, imatorder
            irow(ii:(ii +imatorder -1)) = eltvar(eltptr(inelt):eltptr(inelt+1) -1)
            jcol(ii:(ii +imatorder -1)) = eltvar(eltptr(inelt) +i -1)
            ii = ii +imatorder
         enddo
      enddo

      ! where a_elt brings no contribution, rows and columns are zeroed
      !-----------------------------------------------------------------
      where( abs(a_elt) < EPS_R8 )
         irow = 0
         jcol = 0
      endwhere

      ! the triplet irow, jcol and a_elt is sorted according jcol
      !-----------------------------------------------------------
      call sort_int_1int_1real(g=1, d=ntot, itabref=jcol(1:ntot), itab1=irow(1:ntot), rtab2=a_elt(1:ntot))

      ! column pointer determination for each new value
      !----------------------------------------------
      if (allocated(mat%jptr)) deallocate(mat%jptr) ; allocate( mat%jptr(n +1) ) ; jptr => mat%jptr
      ii = 1
      do
         if (jcol(ii) > 0) exit  ! zeroed terms are ignored
         ii = ii +1
      enddo
      jptr(1) = ii

      do i = 1, n -1
         itmp = jcol(ii)
         do
            ii = ii +1
            if (jcol(ii) /= itmp) exit
         enddo
         jptr(i +1) = ii
      enddo
      jptr(n +1) = ntot +1

      ! columns are already sorted; rows are now sorted for each row value
      !--------------------------------------------------------------------
      do i = 1, n
         i1 = jptr(i)
         i2 = jptr(i +1) -1
         call sort_int_1real(g=1, d=i2 -i1 +1, itabref=irow(i1:i2), rtab1=a_elt(i1:i2))
      enddo

      ! assembly starting from the jcol, irow and a_elt top
      ! for identical matrix locations, a_elts are added
      !-----------------------------------------------------
      innz = 1
      jj  = jptr(1)  ! first non-zero element
      jr1 = jcol(jj)
      ir1 = irow(jj)
      jcol(innz)  = jr1
      irow(innz)  = ir1
      a_elt(innz) = a_elt(jj)

      do j = jj +1, ntot
         jr2 = jcol(j)
         ir2 = irow(j)
         if ( (jr2 /= jr1).or.(ir2 /= ir1) ) then  ! if row or column index has changed
            innz = innz +1                         ! a non-zero term is added
            jcol(innz) = jr2
            irow(innz) = ir2
            a_elt(innz) = a_elt(j)
         else
            a_elt(innz) = a_elt(innz) + a_elt(j)   ! row and column indexes are the same, stiffness terms are added
         endif
         jr1 = jr2                                 ! stores (i-1) and (j-1) for further comparison
         ir1 = ir2
      enddo
      nz = innz
      jcol( nz +1:ntot) = -1
      irow( nz +1:ntot) = -1
      a_elt(nz +1:ntot) = huge(1._R8)

      ! col pointer update
      !----------------------------------------------
      jj      = 1
      jptr(1) = 1
      do j = 1, n -1
         itmp = jcol(jj)
         do
            jj = jj +1
            if (jcol(jj) /= itmp) exit
         enddo
         jptr(j +1) = jj
      enddo
      jptr(n +1) = nz +1

      nullify(eltptr, eltvar, a_elt, irow, jcol, jptr)

   return
   endsubroutine from_elemental_to_assembled

endmodule solver
