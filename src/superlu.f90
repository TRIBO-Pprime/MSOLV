!< author: Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **A SuperLU wrapper**
!  </span>
module sulu_wrapper
use iso_c_binding, only : C_INT, C_FLOAT, C_DOUBLE, C_CHAR, C_PTR, C_F_POINTER, C_NULL_CHAR, C_NULL_PTR
implicit none

private

integer(kind=4), parameter :: I4=4, R4=4, R8=8

!-------------------------------------------------------------------------------
! -------------------   ENUMS --------------------------------------------------
!-------------------------------------------------------------------------------

enum, bind(c) ! trans_t
   enumerator :: NOTRANS
   enumerator :: TRANS
   enumerator :: CONJ
endenum

enum, bind(c) ! fact_t
   enumerator :: DOFACT
   enumerator :: SAMEPATTERN
   enumerator :: SAMEPATTERN_SAMEROWPERM
   enumerator :: FACTORED
endenum

enum, bind(c) ! Stype_t
   enumerator :: SLU_NC     !! *column-wise, no supernode*
   enumerator :: SLU_NCP    !! *column-wise, column-permuted, no supernode*
   enumerator :: SLU_NR     !! *row-wize, no supernode*
   enumerator :: SLU_SC     !! *column-wise, supernode*
   enumerator :: SLU_SCP    !! *supernode, column-wise, permuted*
   enumerator :: SLU_SR     !! *row-wise, supernode*
   enumerator :: SLU_DN     !! *Fortran style column-wise storage for dense matrix*
   enumerator :: SLU_NR_loc !! *distributed compressed row format*
endenum

enum, bind(c) ! Dtype_t
   enumerator :: SLU_S  !! *single*
   enumerator :: SLU_D  !! *double*
   enumerator :: SLU_C  !! *single complex*
   enumerator :: SLU_Z  !! *double complex*
endenum

enum, bind(c) ! Mtype_t
   enumerator :: SLU_GE   !! *general*
   enumerator :: SLU_TRLU !! *lower triangular, unit diagonal*
   enumerator :: SLU_TRUU !! *upper triangular, unit diagonal*
   enumerator :: SLU_TRL  !! *lower triangular*
   enumerator :: SLU_TRU  !! *upper triangular*
   enumerator :: SLU_SYL  !! *symmetric, store lower half*
   enumerator :: SLU_SYU  !! *symmetric, store upper half*
   enumerator :: SLU_HEL  !! *Hermitian, store lower half*
   enumerator :: SLU_HEU  !! *Hermitian, store upper half*
endenum

!-------------------------------------------------------------------------------
! -------------------   DERIVED TYPES   ----------------------------------------
!-------------------------------------------------------------------------------

type, bind(c) :: LU_STACK_T
   integer(kind=C_INT) :: size
   integer(kind=C_INT) :: used
   integer(kind=C_INT) :: top1
   integer(kind=C_INT) :: top2
   type(C_PTR)         :: array
endtype LU_STACK_T

type, bind(c) :: EXPHEADER
   integer(kind=C_INT) :: size
   type(C_PTR)    :: mem
endtype EXPHEADER

type, bind(c) :: GLOBALLU_T
   integer(kind=C_INT) :: xsup
   integer(kind=C_INT) :: supno
   integer(kind=C_INT) :: lsub
   integer(kind=C_INT) :: xlsub
   type(C_PTR)         :: lusup
   integer(kind=C_INT) :: xlusup
   type(C_PTR)         :: ucol
   integer(kind=C_INT) :: usub
   integer(kind=C_INT) :: xusub
   integer(kind=C_INT) :: nzlmax
   integer(kind=C_INT) :: nzumax
   integer(kind=C_INT) :: nzlumax
   integer(kind=C_INT) :: MemModel=0
   integer(kind=C_INT) :: num_expansions
   type(EXPHEADER)     :: expanders
   type(LU_STACK_T)    :: stack
endtype GLOBALLU_T

type, bind(c) :: SUPERLUSTAT_T
   type(C_PTR)         :: panel_histo    !! *histogram of panel size distribution*
   type(C_PTR)         :: utime          !! *running time at various phases*
   type(C_PTR)         :: ops            !! *operation count at various phases*
   integer(kind=C_INT) :: TinyPivots     !! *number of tiny pivots*
   integer(kind=C_INT) :: RefineSteps    !! *number of iterative refinement steps*
   integer(kind=C_INT) :: expansions     !! *number of memory expansions*
endtype SUPERLUSTAT_T

type, bind(c) :: MEM_USAGE_T
   real(kind=C_FLOAT) :: for_lu
   real(kind=C_FLOAT) :: total_needed
endtype MEM_USAGE_T

!< @note
! Stype == ```SLU_NC``` (Also known as Harwell-Boeing sparse matrix format)
!
! Zero-based indexing is used; colptr[] has ncol+1 entries, the last one pointing beyond the last column,
! so that colptr[ncol] = nnz.
! @endnote
type, bind(c) :: NCFORMAT 
   integer(kind=C_INT) :: nnz    !! *number of nonzeros in the matrix*
   type(C_PTR)         :: nzval  !! *pointer to array of nonzero values, packed by column*
   type(C_PTR)         :: rowind !! *pointer to array of row indices of the nonzeros*
   type(C_PTR)         :: colptr !! *pointer to array of beginning of columns in nzval[] and rowind[]*
endtype NCFORMAT

type, bind(c) :: SUPERMATRIX
   integer(kind=C_INT) :: Stype    !! *Storage type: interprets the storage structure pointed to by Store*
   integer(kind=C_INT) :: Dtype    !! *Data type*
   integer(kind=C_INT) :: Mtype    !! *Matrix type: describes the mathematical property of the matrix*
   integer(kind=C_INT) :: nrow     !! *number of rows*
   integer(kind=C_INT) :: ncol     !! *number of columns*
   type(C_PTR)         :: Store    !! *pointer to the actual storage of the matrix, here, pointer to [[NCformat]]*
endtype SUPERMATRIX

type, bind(c) :: SUPERLU_OPTIONS_T
   !/*
   ! *-- This contains the options used to control the solution process.
   ! *
   ! * Fact   (fact_t)
   ! *        Specifies whether or not the factored form of the matrix
   ! *        A is supplied on entry, and if not, how the matrix A should
   ! *        be factorizaed.
   ! *        = DOFACT: The matrix A will be factorized from scratch, and the
   ! *             factors will be stored in L and U.
   ! *        = SamePattern: The matrix A will be factorized assuming
   ! *             that a factorization of a matrix with the same sparsity
   ! *             pattern was performed prior to this one. Therefore, this
   ! *             factorization will reuse column permutation vector
   ! *             ScalePermstruct->perm_c and the column elimination tree
   ! *             LUstruct->etree.
   ! *        = SamePattern_SameRowPerm: The matrix A will be factorized
   ! *             assuming that a factorization of a matrix with the same
   ! *             sparsity   pattern and similar numerical values was performed
   ! *             prior to this one. Therefore, this factorization will reuse
   ! *             both row and column scaling factors R and C, both row and
   ! *             column permutation vectors perm_r and perm_c, and the
   ! *             L & U data structures set up from the previous factorization.
   !! *        = FACTORED: On entry, L, U, perm_r and perm_c contain the
   ! *              factored form of A. If DiagScale is not NOEQUIL, the matrix
   ! *              A has been equilibrated with scaling factors R and C.
   ! *
   ! * Equil  (yes_no_t)
   ! *        Specifies whether to equilibrate the system (scale A's row and
   ! *        columns to have unit norm).
   ! *
   ! * ColPerm (colperm_t)
   ! *        Specifies what type of column permutation to use to reduce fill.
   ! *        = NATURAL: use the natural ordering
   ! *        = MMD_ATA: use minimum degree ordering on structure of A'*A
   ! *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
   ! *        = COLAMD: use approximate minimum degree column ordering
   ! *        = MY_PERMC: use the ordering specified by the user
   ! *
   ! * Trans  (trans_t)
   ! *        Specifies the form of the system of equations:
   ! *        = NOTRANS: A * X = B        (No transpose)
   ! *        = TRANS:   A**T * X = B     (Transpose)
   ! *        = CONJ:    A**H * X = B     (Transpose)
   ! *
   ! * IterRefine (IterRefine_t)
   ! *        Specifies whether to perform iterative refinement.
   ! *        = NO: no iterative refinement
   ! *        = SLU_SINGLE: perform iterative refinement in single precision
   ! *        = SLU_DOUBLE: perform iterative refinement in double precision
   ! *        = SLU_EXTRA: perform iterative refinement in extra precision
   ! *
   ! * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
   ! *        Specifies the threshold used for a diagonal entry to be an
   ! *        acceptable pivot.
   ! *
   ! * SymmetricMode (yest_no_t)
   ! *        Specifies whether to use symmetric mode. Symmetric mode gives
   ! *        preference to diagonal pivots, and uses an (A'+A)-based column
   ! *        permutation algorithm.
   ! *
   ! * PivotGrowth (yes_no_t)
   ! *        Specifies whether to compute the reciprocal pivot growth.
   ! *
   ! * ConditionNumber (ues_no_t)
   ! *        Specifies whether to compute the reciprocal condition number.
   ! *
   ! * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
   ! *        Specifies whether to permute rows of the original matrix.
   ! *        = NO: not to permute the rows
   ! *        = LargeDiag: make the diagonal large relative to the off-diagonal
   ! *        = MY_PERMR: use the permutation given by the user
   ! *
   ! * ILU_DropRule (int)
   ! *        Specifies the dropping rule:
   ! *     = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
   ! *     = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
   ! *     = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
   ! *                     p = gamma * nnz(A(:,j)).
   ! *     = DROP_AREA:    Variation of ILUTP, for j-th column, use
   ! *                     nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
   ! *     = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
   ! *                     If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
   ! *                         tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
   ! *                     Otherwise
   ! *                         tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
   ! *                     tau_U(j) uses the similar rule.
   ! *                     NOTE: the thresholds used by L and U are separate.
   ! *     = DROP_INTERP:  Compute the second dropping threshold by
   ! *                     interpolation instead of sorting (default).
   ! *                     In this case, the actual fill ratio is not
   ! *                     guaranteed to be smaller than gamma.
   ! *                     Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
   ! *    ( Default: DROP_BASIC | DROP_AREA )
   ! *
   ! * ILU_DropTol (double)
   ! *        numerical threshold for dropping.
   ! *
   ! * ILU_FillFactor (double)
   ! *        Gamma in the secondary dropping.
   ! *
   ! * ILU_Norm (norm_t)
   ! *        Specify which norm to use to measure the row size in a
   ! *        supernode: infinity-norm, 1-norm, or 2-norm.
   ! *
   ! * ILU_FillTol (double)
   ! *        numerical threshold for zero pivot perturbation.
   ! *
   ! * ILU_MILU (milu_t)
   ! *        Specifies which version of MILU to use.
   ! *
   ! * ILU_MILU_Dim (double)
   ! *        Dimension of the PDE if available.
   ! *
   ! * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
   ! *        Specifies whether to replace the tiny diagonals by
   ! *        sqrt(epsilon)*||A|| during LU factorization.
   ! *
   ! * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
   ! *        Specifies whether the initialization has been performed to the
   ! *        triangular solve.
   ! *
   ! * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
   ! *        Specifies whether the initialization has been performed to the
   ! *        sparse matrix-vector multiplication routine needed in iterative
   ! *        refinement.
   ! *
   ! * PrintStat (yes_no_t)
   ! *        Specifies whether to print the solver's statistics.
   ! */
   integer(kind=C_INT) :: Fact
   integer(kind=C_INT) :: Equil
   integer(kind=C_INT) :: ColPerm
   integer(kind=C_INT) :: Trans
   integer(kind=C_INT) :: IterRefine
   real(kind=C_DOUBLE) :: DiagPivotThresh
   integer(kind=C_INT) :: SymmetricMode
   integer(kind=C_INT) :: PivotGrowth
   integer(kind=C_INT) :: ConditionNumber
   integer(kind=C_INT) :: RowPerm
   integer(kind=C_INT) :: ILU_DropRule
   real(kind=C_DOUBLE) :: ILU_DropTol
   real(kind=C_DOUBLE) :: ILU_FillFactor
   integer(kind=C_INT) :: ILU_Norm
   real(kind=C_DOUBLE) :: ILU_FillTol
   integer(kind=C_INT) :: ILU_MILU
   real(kind=C_DOUBLE) :: ILU_MILU_Dim
   integer(kind=C_INT) :: ParSymbFact
   integer(kind=C_INT) :: ReplaceTinyPivot
   integer(kind=C_INT) :: SolveInitialized
   integer(kind=C_INT) :: RefineInitialized
   integer(kind=C_INT) :: PrintStat
   integer(kind=C_INT) :: nnzL, nnzU      !! *used to store nnzs for now*
   integer(kind=C_INT) :: num_lookaheads  !! *num of levels in look-ahead*
   integer(kind=C_INT) :: lookahead_etree !! *use etree computed from the serial symbolic factorization*
   integer(kind=C_INT) :: SymPattern      !! *symmetric factorization*
endtype SUPERLU_OPTIONS_T

!-------------------------------------------------------------------------------
! -------------------   SULU GLOBAL TYPE   -------------------------------------
!-------------------------------------------------------------------------------

type SULU_ENV
!! <span style="color:green">Global type for *SuperLU* which covers all the stuff needed</span>
   integer(kind=C_INT) :: n      !! *system size*
   integer(kind=C_INT) :: nrhs   !! *number of right hand sides*
   integer(kind=C_INT) :: nz     !! *number on non-zero entries*
   integer(kind=C_INT) :: info   !! *info returned by [[dgssvx]]*
   integer(kind=C_INT) :: lwork  !! *size of workspace, not used here*

   logical(kind=I4)    :: first  !! *if ```false``` the system has been factorized at least once*

   real(kind=R8), dimension(:), pointer     :: b   !! *right hand side: points to [[MAT_SOLV:b]]*
   real(kind=R8), allocatable, dimension(:) :: x   !! *solution*

   real(kind=R8),    dimension(:), pointer :: a_elt   !! *CC system matrix: points to [[MAT_SOLV:a_elt]]*
   integer(kind=I4), dimension(:), pointer :: irow    !! *matrix line of an a_elt element: points to [[MAT_SOLV:irow]]*
   integer(kind=I4), dimension(:), pointer :: jptr    !! *matrix column pointers: points to [[MAT_SOLV:jptr]]*

   real(kind=C_DOUBLE), allocatable, dimension(:) :: ferr   !! *estimated forward error bound for each solution vector*
   real(kind=C_DOUBLE), allocatable, dimension(:) :: berr   !! *componentwise relative backward error of each solution*
   real(kind=C_DOUBLE), allocatable, dimension(:) :: RR     !! *row scale factors for A *
   real(kind=C_DOUBLE), allocatable, dimension(:) :: CC     !!*column scale factors for A*
   real(kind=C_DOUBLE), allocatable, dimension(:) :: rpg    !! *reciprocal pivot growth factor*
   real(kind=C_DOUBLE), allocatable, dimension(:) :: rcond  !!*estimate of the reciprocal condition number of the matrix A*

   integer(kind=C_INT), allocatable, dimension(:) :: perm_c !!*If A->Stype = ```SLU_NC```, Column permutation vector of size A->ncol*
   integer(kind=C_INT), allocatable, dimension(:) :: perm_r !!*If A->Stype = ```SLU_NC```, row permutation vector of size A->nrow*
   integer(kind=C_INT), allocatable, dimension(:) :: etree  !! *Elimination tree*

   character(len=1, kind=C_CHAR) :: equed !! *form of equilibration*

   type(C_PTR) :: work                    !! *User supplied workspace*

   type(SUPERLU_OPTIONS_T) :: options     !! *LU controls*

   type(SUPERMATRIX) :: sma !! *Matrix A in A*X=B*
   type(SUPERMATRIX) :: smb !! *On entry, the right hand side matrix*
   type(SUPERMATRIX) :: smx !! *olution matrix to the original system*
   type(SUPERMATRIX) :: sml !! *factor L from the factorization*
   type(SUPERMATRIX) :: smu !! *factor U from the factorization*

   type(SUPERLUSTAT_T) :: stat !! *statistics on runtime and floating-point operation count*

   type(GLOBALLU_T)    :: Glu  !! *first, an output with the whole stuff LU; next, an input for other resolutions with same sparsity*

   type(MEM_USAGE_T)   :: mem_usage !! *memory usage statistics*
endtype SULU_ENV

!-------------------------------------------------------------------------------
! -------------------   INTERFACES   -------------------------------------------
!-------------------------------------------------------------------------------

interface

   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine Destroy_SuperNode_Matrix(A) & !
   bind(c, name="Destroy_SuperNode_Matrix")
   import :: SUPERMATRIX
   implicit none
   type(SUPERMATRIX), intent(in) :: A
   endsubroutine Destroy_SuperNode_Matrix
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine Destroy_SuperMatrix_Store(A) & !
   bind(c, name="Destroy_SuperMatrix_Store")
   import :: SUPERMATRIX
   implicit none
   type(SUPERMATRIX), intent(in) :: A
   endsubroutine Destroy_SuperMatrix_Store
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine Destroy_CompCol_Matrix(A) & !
   bind(c, name="Destroy_CompCol_Matrix")
   import :: SUPERMATRIX
   implicit none
   type(SUPERMATRIX), intent(in) :: A
   endsubroutine Destroy_CompCol_Matrix
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine Destroy_Dense_Matrix(A) & !
   bind(c, name="Destroy_Dense_Matrix")
   import :: SUPERMATRIX
   implicit none
   type(SUPERMATRIX), intent(in) :: A
   endsubroutine Destroy_Dense_Matrix
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine StatInit(stat) & !
   bind(c, name="StatInit")
   import :: SUPERLUSTAT_T
   implicit none
   type(SUPERLUSTAT_T), intent(out) :: stat
   endsubroutine StatInit
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine StatFree(stat) & !
   bind(c, name="StatFree")
   import :: SUPERLUSTAT_T
   implicit none
   type(SUPERLUSTAT_T), intent(in) :: stat
   endsubroutine StatFree
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine dCreate_CompCol_Matrix(  A,       & ! out SuperMatrix
                                       m,       & ! in int
                                       n,       & ! in int
                                       nnz,     & ! in int
                                       nzval,   & ! in double dimension()
                                       rowind,  & ! in int dimension()
                                       colptr,  & ! in int dimension()
                                       stype,   & ! in int
                                       dtype,   & ! in int
                                       mtype    & ! in int
                                    ) &
   bind(c, name="dCreate_CompCol_Matrix")
   use, intrinsic :: iso_c_binding, only : C_INT, C_DOUBLE
   import         :: SUPERMATRIX
   implicit none
   type(SUPERMATRIX),             intent(out) :: A
   integer(kind=C_INT),    value, intent(in)  :: m
   integer(kind=C_INT),    value, intent(in)  :: n
   integer(kind=C_INT),    value, intent(in)  :: nnz
   real(kind=C_DOUBLE),           intent(in)  :: nzval(*)
   integer(kind=C_INT),           intent(in)  :: rowind(*)
   integer(kind=C_INT),           intent(in)  :: colptr(*)
   integer(kind=C_INT),    value, intent(in)  :: stype
   integer(kind=C_INT),    value, intent(in)  :: dtype
   integer(kind=C_INT),    value, intent(in)  :: mtype
   endsubroutine dCreate_CompCol_Matrix
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine dCreate_Dense_Matrix( BX,      & ! out SuperMatrix
                                    m,       & ! in int
                                    n,       & ! in int
                                    x,       & ! in double dimension()
                                    ldx,     & ! in int
                                    stype,   & ! in int
                                    dtype,   & ! in int
                                    mtype    & ! in int
                                    ) &
   bind(c, name='dCreate_Dense_Matrix')
   use, intrinsic :: iso_c_binding, only : C_INT, C_DOUBLE
   import SUPERMATRIX
   implicit none
   type(SUPERMATRIX),              intent(out) :: BX
   integer(kind=C_INT),    value,  intent(in)  :: m
   integer(kind=C_INT),    value,  intent(in)  :: n
   real(kind=C_DOUBLE),            intent(in)  :: x(*)
   integer(kind=C_INT),    value,  intent(in)  :: ldx
   integer(kind=C_INT),    value,  intent(in)  :: stype
   integer(kind=C_INT),    value,  intent(in)  :: dtype
   integer(kind=C_INT),    value,  intent(in)  :: mtype
   endsubroutine dCreate_Dense_Matrix
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine set_default_options(options) & !
   bind(c, name='set_default_options')
   use, intrinsic :: iso_c_binding
   import SUPERLU_OPTIONS_T
   implicit none
   type(SUPERLU_OPTIONS_T), intent(inout) :: options
   endsubroutine set_default_options
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   !/*! Arguments
   ! *
   ! * <pre>
   ! * Purpose
   ! * =======
   ! *
   ! * DGSSVX solves the system of linear equations A*X=B or A'*X=B, using
   ! * the LU factorization from dgstrf(). Error bounds on the solution and
   ! * a condition estimate are also provided. It performs the following steps:
   ! *
   ! *   1. If A is stored column-wise (A->Stype = SLU_NC):
   ! *
   ! *      1.1. If options->Equil = YES, scaling factors are computed to
   ! *           equilibrate the system:
   ! *           options->Trans = NOTRANS:
   ! *               diag(R)*A*diag(C) *inv(diag(C))*X = diag(R)*B
   ! *           options->Trans = TRANS:
   ! *               (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
   ! *           options->Trans = CONJ:
   ! *               (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
   ! *           Whether or not the system will be equilibrated depends on the
   ! *           scaling of the matrix A, but if equilibration is used, A is
   ! *           overwritten by diag(R)*A*diag(C) and B by diag(R)*B
   ! *           (if options->Trans=NOTRANS) or diag(C)*B (if options->Trans
   ! *           = TRANS or CONJ).
   ! *
   ! *      1.2. Permute columns of A, forming A*Pc, where Pc is a permutation
   ! *           matrix that usually preserves sparsity.
   ! *           For more details of this step, see sp_preorder.c.
   ! *
   ! *      1.3. If options->Fact != FACTORED, the LU decomposition is used to
   ! *           factor the matrix A (after equilibration if options->Equil = YES)
   ! *           as Pr*A*Pc = L*U, with Pr determined by partial pivoting.
   ! *
   ! *      1.4. Compute the reciprocal pivot growth factor.
   ! *
   ! *      1.5. If some U(i,i) = 0, so that U is exactly singular, then the
   ! *           routine returns with info = i. Otherwise, the factored form of
   ! *           A is used to estimate the condition number of the matrix A. If
   ! *           the reciprocal of the condition number is less than machine
   ! *           precision, info = A->ncol+1 is returned as a warning, but the
   ! *           routine still goes on to solve for X and computes error bounds
   ! *           as described below.
   ! *
   ! *      1.6. The system of equations is solved for X using the factored form
   ! *           of A.
   ! *
   ! *      1.7. If options->IterRefine != NOREFINE, iterative refinement is
   ! *           applied to improve the computed solution matrix and calculate
   ! *           error bounds and backward error estimates for it.
   ! *
   ! *      1.8. If equilibration was used, the matrix X is premultiplied by
   ! *           diag(C) (if options->Trans = NOTRANS) or diag(R)
   ! *           (if options->Trans = TRANS or CONJ) so that it solves the
   ! *           original system before equilibration.
   ! *
   ! *   2. If A is stored row-wise (A->Stype = SLU_NR), apply the above algorithm
   ! *      to the transpose of A:
   ! *
   ! *      2.1. If options->Equil = YES, scaling factors are computed to
   ! *           equilibrate the system:
   ! *           options->Trans = NOTRANS:
   ! *               diag(R)*A*diag(C) *inv(diag(C))*X = diag(R)*B
   ! *           options->Trans = TRANS:
   ! *               (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
   ! *           options->Trans = CONJ:
   ! *               (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
   ! *           Whether or not the system will be equilibrated depends on the
   ! *           scaling of the matrix A, but if equilibration is used, A' is
   ! *           overwritten by diag(R)*A'*diag(C) and B by diag(R)*B
   ! *           (if trans='N') or diag(C)*B (if trans = 'T' or 'C').
   ! *
   ! *      2.2. Permute columns of transpose(A) (rows of A),
   ! *           forming transpose(A)*Pc, where Pc is a permutation matrix that
   ! *           usually preserves sparsity.
   ! *           For more details of this step, see sp_preorder.c.
   ! *
   ! *      2.3. If options->Fact != FACTORED, the LU decomposition is used to
   ! *           factor the transpose(A) (after equilibration if
   ! *           options->Fact = YES) as Pr*transpose(A)*Pc = L*U with the
   ! *           permutation Pr determined by partial pivoting.
   ! *
   ! *      2.4. Compute the reciprocal pivot growth factor.
   ! *
   ! *      2.5. If some U(i,i) = 0, so that U is exactly singular, then the
   ! *           routine returns with info = i. Otherwise, the factored form
   ! *           of transpose(A) is used to estimate the condition number of the
   ! *           matrix A. If the reciprocal of the condition number
   ! *           is less than machine precision, info = A->nrow+1 is returned as
   ! *           a warning, but the routine still goes on to solve for X and
   ! *           computes error bounds as described below.
   ! *
   ! *      2.6. The system of equations is solved for X using the factored form
   ! *           of transpose(A).
   ! *
   ! *      2.7. If options->IterRefine != NOREFINE, iterative refinement is
   ! *           applied to improve the computed solution matrix and calculate
   ! *           error bounds and backward error estimates for it.
   ! *
   ! *      2.8. If equilibration was used, the matrix X is premultiplied by
   ! *           diag(C) (if options->Trans = NOTRANS) or diag(R)
   ! *           (if options->Trans = TRANS or CONJ) so that it solves the
   ! *           original system before equilibration.
   ! *
   ! *   See supermatrix.h for the definition of 'SuperMatrix' structure.
   ! *
   ! * Arguments
   ! * =========
   ! *
   ! * options (input) superlu_options_t*
   ! *         The structure defines the input parameters to control
   ! *         how the LU decomposition will be performed and how the
   ! *         system will be solved.
   ! *
   ! * A       (input/output) SuperMatrix*
   ! *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
   ! *         of the linear equations is A->nrow. Currently, the type of A can be:
   ! *         Stype = SLU_NC or SLU_NR, Dtype = SLU_D, Mtype = SLU_GE.
   ! *         In the future, more general A may be handled.
   ! *
   ! *         On entry, If options->Fact = FACTORED and equed is not 'N',
   ! *         then A must have been equilibrated by the scaling factors in
   ! *         R and/or C.
   ! *         On exit, A is not modified if options->Equil = NO, or if
   ! *         options->Equil = YES but equed = 'N' on exit.
   ! *         Otherwise, if options->Equil = YES and equed is not 'N',
   ! *         A is scaled as follows:
   ! *         If A->Stype = SLU_NC:
   ! *           equed = 'R':  A := diag(R) * A
   ! *           equed = 'C':  A := A * diag(C)
   ! *           equed = 'B':  A := diag(R) * A * diag(C).
   ! *         If A->Stype = SLU_NR:
   ! *           equed = 'R':  transpose(A) := diag(R) * transpose(A)
   ! *           equed = 'C':  transpose(A) := transpose(A) * diag(C)
   ! *           equed = 'B':  transpose(A) := diag(R) * transpose(A) * diag(C).
   ! *
   ! * perm_c  (input/output) int*
   ! *         If A->Stype = SLU_NC, Column permutation vector of size A->ncol,
   ! *         which defines the permutation matrix Pc; perm_c[i] = j means
   ! *         column i of A is in position j in A*Pc.
   ! *         On exit, perm_c may be overwritten by the product of the input
   ! *         perm_c and a permutation that postorders the elimination tree
   ! *         of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
   ! *         is already in postorder.
   ! *
   ! *         If A->Stype = SLU_NR, column permutation vector of size A->nrow,
   ! *         which describes permutation of columns of transpose(A)
   ! *         (rows of A) as described above.
   ! *
   ! * perm_r  (input/output) int*
   ! *         If A->Stype = SLU_NC, row permutation vector of size A->nrow,
   ! *         which defines the permutation matrix Pr, and is determined
   ! *         by partial pivoting.  perm_r[i] = j means row i of A is in
   ! *         position j in Pr*A.
   ! *
   ! *         If A->Stype = SLU_NR, permutation vector of size A->ncol, which
   ! *         determines permutation of rows of transpose(A)
   ! *         (columns of A) as described above.
   ! *
   ! *         If options->Fact = SamePattern_SameRowPerm, the pivoting routine
   ! *         will try to use the input perm_r, unless a certain threshold
   ! *         criterion is violated. In that case, perm_r is overwritten by a
   ! *         new permutation determined by partial pivoting or diagonal
   ! *         threshold pivoting.
   ! *         Otherwise, perm_r is output argument.
   ! *
   ! * etree   (input/output) int*,  dimension (A->ncol)
   ! *         Elimination tree of Pc'*A'*A*Pc.
   ! *         If options->Fact != FACTORED and options->Fact != DOFACT,
   ! *         etree is an input argument, otherwise it is an output argument.
   ! *         Note: etree is a vector of parent pointers for a forest whose
   ! *         vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
   ! *
   ! * equed   (input/output) char*
   ! *         Specifies the form of equilibration that was done.
   ! *         = 'N': No equilibration.
   ! *         = 'R': Row equilibration, i.e., A was premultiplied by diag(R).
   ! *         = 'C': Column equilibration, i.e., A was postmultiplied by diag(C).
   ! *         = 'B': Both row and column equilibration, i.e., A was replaced
   ! *                by diag(R)*A*diag(C).
   ! *         If options->Fact = FACTORED, equed is an input argument,
   ! *         otherwise it is an output argument.
   ! *
   ! * R       (input/output) double*, dimension (A->nrow)
   ! *         The row scale factors for A or transpose(A).
   ! *         If equed = 'R' or 'B', A (if A->Stype = SLU_NC) or transpose(A)
   ! *             (if A->Stype = SLU_NR) is multiplied on the left by diag(R).
   ! *         If equed = 'N' or 'C', R is not accessed.
   ! *         If options->Fact = FACTORED, R is an input argument,
   ! *             otherwise, R is output.
   ! *         If options->zFact = FACTORED and equed = 'R' or 'B', each element
   ! *             of R must be positive.
   ! *
   ! * C       (input/output) double*, dimension (A->ncol)
   ! *         The column scale factors for A or transpose(A).
   ! *         If equed = 'C' or 'B', A (if A->Stype = SLU_NC) or transpose(A)
   ! *             (if A->Stype = SLU_NR) is multiplied on the right by diag(C).
   ! *         If equed = 'N' or 'R', C is not accessed.
   ! *         If options->Fact = FACTORED, C is an input argument,
   ! *             otherwise, C is output.
   ! *         If options->Fact = FACTORED and equed = 'C' or 'B', each element
   ! *             of C must be positive.
   ! *
   ! * L       (output) SuperMatrix*
   ! *         The factor L from the factorization
   ! *             Pr*A*Pc=L*U              (if A->Stype SLU_= NC) or
   ! *             Pr*transpose(A)*Pc=L*U   (if A->Stype = SLU_NR).
   ! *         Uses compressed row subscripts storage for supernodes, i.e.,
   ! *         L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
   ! *
   ! * U       (output) SuperMatrix*
   ! *         The factor U from the factorization
   ! *             Pr*A*Pc=L*U              (if A->Stype = SLU_NC) or
   ! *             Pr*transpose(A)*Pc=L*U   (if A->Stype = SLU_NR).
   ! *         Uses column-wise storage scheme, i.e., U has types:
   ! *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
   ! *
   ! * work    (workspace/output) void*, size (lwork) (in bytes)
   ! *         User supplied workspace, should be large enough
   ! *         to hold data structures for factors L and U.
   ! *         On exit, if fact is not 'F', L and U point to this array.
   ! *
   ! * lwork   (input) int
   ! *         Specifies the size of work array in bytes.
   ! *         = 0:  allocate space internally by system malloc;
   ! *         > 0:  use user-supplied work array of length lwork in bytes,
   ! *               returns error if space runs out.
   ! *         = -1: the routine guesses the amount of space needed without
   ! *               performing the factorization, and returns it in
   ! *               mem_usage->total_needed; no other side effects.
   ! *
   ! *         See argument 'mem_usage' for memory usage statistics.
   ! *
   ! * B       (input/output) SuperMatrix*
   ! *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
   ! *         On entry, the right hand side matrix.
   ! *         If B->ncol = 0, only LU decomposition is performed, the triangular
   ! *                         solve is skipped.
   ! *         On exit,
   ! *            if equed = 'N', B is not modified; otherwise
   ! *            if A->Stype = SLU_NC:
   ! *               if options->Trans = NOTRANS and equed = 'R' or 'B',
   ! *                  B is overwritten by diag(R)*B;
   ! *               if options->Trans = TRANS or CONJ and equed = 'C' of 'B',
   ! *                  B is overwritten by diag(C)*B;
   ! *            if A->Stype = SLU_NR:
   ! *               if options->Trans = NOTRANS and equed = 'C' or 'B',
   ! *                  B is overwritten by diag(C)*B;
   ! *               if options->Trans = TRANS or CONJ and equed = 'R' of 'B',
   ! *                  B is overwritten by diag(R)*B.
   ! *
   ! * X       (output) SuperMatrix*
   ! *         X has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
   ! *         If info = 0 or info = A->ncol+1, X contains the solution matrix
   ! *         to the original system of equations. Note that A and B are modified
   ! *         on exit if equed is not 'N', and the solution to the equilibrated
   ! *         system is inv(diag(C))*X if options->Trans = NOTRANS and
   ! *         equed = 'C' or 'B', or inv(diag(R))*X if options->Trans = 'T' or 'C'
   ! *         and equed = 'R' or 'B'.
   ! *
   ! * recip_pivot_growth (output) double*
   ! *         The reciprocal pivot growth factor max_j( norm(A_j)/norm(U_j) ).
   ! *         The infinity norm is used. If recip_pivot_growth is much less
   ! *         than 1, the stability of the LU factorization could be poor.
   ! *
   ! * rcond   (output) double*
   ! *         The estimate of the reciprocal condition number of the matrix A
   ! *         after equilibration (if done). If rcond is less than the machine
   ! *         precision (in particular, if rcond = 0), the matrix is singular
   ! *         to working precision. This condition is indicated by a return
   ! *         code of info > 0.
   ! *
   ! * FERR    (output) double*, dimension (B->ncol)
   ! *         The estimated forward error bound for each solution vector
   ! *         X(j) (the j-th column of the solution matrix X).
   ! *         If XTRUE is the true solution corresponding to X(j), FERR(j)
   ! *         is an estimated upper bound for the magnitude of the largest
   ! *         element in (X(j) - XTRUE) divided by the magnitude of the
   ! *         largest element in X(j).  The estimate is as reliable as
   ! *         the estimate for RCOND, and is almost always a slight
   ! *         overestimate of the true error.
   ! *         If options->IterRefine = NOREFINE, ferr = 1.0.
   ! *
   ! * BERR    (output) double*, dimension (B->ncol)
   ! *         The componentwise relative backward error of each solution
   ! *         vector X(j) (i.e., the smallest relative change in
   ! *         any element of A or B that makes X(j) an exact solution).
   ! *         If options->IterRefine = NOREFINE, berr = 1.0.
   ! *
   ! * Glu      (input/output) GlobalLU_t *
   ! *          If options->Fact == SamePattern_SameRowPerm, it is an input;
   ! *              The matrix A will be factorized assuming that a
   ! *              factorization of a matrix with the same sparsity pattern
   ! *              and similar numerical values was performed prior to this one.
   ! *              Therefore, this factorization will reuse both row and column
   ! *              scaling factors R and C, both row and column permutation
   ! *              vectors perm_r and perm_c, and the L & U data structures
   ! *              set up from the previous factorization.
   ! *          Otherwise, it is an output.
   ! *
   ! * mem_usage (output) mem_usage_t*
   ! *         Record the memory usage statistics, consisting of following fields:
   ! *         - for_lu (float)
   ! *           The amount of space used in bytes for L\U data structures.
   ! *         - total_needed (float)
   ! *           The amount of space needed in bytes to perform factorization.
   ! *         - expansions (int)
   ! *           The number of memory expansions during the LU factorization.
   ! *
   ! * stat   (output) SuperLUStat_t*
   ! *        Record the statistics on runtime and floating-point operation count.
   ! *        See slu_util.h for the definition of 'SuperLUStat_t'.
   ! *
   ! * info    (output) int*
   ! *         = 0: successful exit
   ! *         < 0: if info = -i, the i-th argument had an illegal value
   ! *         > 0: if info = i, and i is
   ! *              <= A->ncol: U(i,i) is exactly zero. The factorization has
   ! *                    been completed, but the factor U is exactly
   ! *                    singular, so the solution and error bounds
   ! *                    could not be computed.
   ! *              = A->ncol+1: U is nonsingular, but RCOND is less than machine
   ! *                    precision, meaning that the matrix is singular to
   ! *                    working precision. Nevertheless, the solution and
   ! *                    error bounds are computed because there are a number
   ! *                    of situations where the computed solution can be more
   ! *                    accurate than the value of RCOND would suggest.
   ! *              > A->ncol+1: number of bytes allocated when memory allocation
   ! *                    failure occurred, plus A->ncol.
   ! * </pre>
   ! */
   subroutine dgssvx(& ! argument                 |     type           |   C def           |  C call
                        options,                & ! superlu_options_t   *options             &options
                        A,                      & ! SuperMatrix         *A                   &A1
                        perm_c,                 & ! int                 *perm_c               perm_c
                        perm_r,                 & ! int                 *perm_r               perm_r
                        etree,                  & ! int                 *etree                etree
                        equed,                  & ! char                *equed                equed
                        R,                      & ! double              *R                    R
                        C,                      & ! double              *C                    C
                        L,                      & ! SuperMatrix         *L                   &L
                        U,                      & ! SuperMatrix         *U                   &U
                        work,                   & ! void                *work                 work
                        lwork,                  & ! int                  lwork                lwork
                        B,                      & ! SuperMatrix         *B                   &B1
                        X,                      & ! SuperMatrix         *X                   &X
                        recip_pivot_growth,     & ! double              *recip_pivot_growth  &rpg
                        rcond,                  & ! double              *rcond               &rcond
                        ferr,                   & ! double              *ferr                 ferr
                        berr,                   & ! double              *berr                 berr
                        Glu,                    & ! GlobalLU_t          *Glu                 &Glu
                        mem_usage,              & ! mem_usage_t         *mem_usage           &smem_usage
                        stat,                   & ! SuperLUStat_t       *stat                &stat
                        info                    & ! int                 *info                &info
                     ) &
   bind(c, name='dgssvx')
   use, intrinsic :: iso_c_binding, only : C_INT, C_CHAR, C_DOUBLE, C_PTR
   import :: SUPERLU_OPTIONS_T, SUPERMATRIX, SUPERLUSTAT_T, MEM_USAGE_T, GLOBALLU_T
   implicit none
   type(SUPERLU_OPTIONS_T),        intent(in)     :: options
   type(SUPERMATRIX),              intent(inout)  :: A
   integer(kind=C_INT),            intent(inout)  :: perm_c(*)
   integer(kind=C_INT),            intent(inout)  :: perm_r(*)
   integer(kind=C_INT),            intent(inout)  :: etree(*)
   character(kind=C_CHAR),         intent(inout)  :: equed(*)
   real(kind=C_DOUBLE),            intent(inout)  :: R(*)
   real(kind=C_DOUBLE),            intent(inout)  :: C(*)
   type(SUPERMATRIX),              intent(inout)  :: L
   type(SUPERMATRIX),              intent(inout)  :: U
   type(C_PTR),                    intent(out)    :: work
   integer(kind=C_INT),    value,  intent(in)     :: lwork
   type(SUPERMATRIX),              intent(inout)  :: B
   type(SUPERMATRIX),              intent(out)    :: X
   real(kind=C_DOUBLE),            intent(out)    :: recip_pivot_growth(*)
   real(kind=C_DOUBLE),            intent(out)    :: rcond(*)
   real(kind=C_DOUBLE),            intent(out)    :: ferr(*)
   real(kind=C_DOUBLE),            intent(out)    :: berr(*)
   type(GLOBALLU_T),               intent(inout)  :: Glu
   type(MEM_USAGE_T),              intent(out)    :: mem_usage
   type(SUPERLUSTAT_T),            intent(out)    :: stat
   integer(kind=C_INT),            intent(out)    :: info
   endsubroutine dgssvx
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   subroutine StatPrint(stat) & !
   bind(c, name="StatPrint")
   import :: SUPERLUSTAT_T
   implicit none
   type(SUPERLUSTAT_T), intent(in) :: stat
   endsubroutine StatPrint
   !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

endinterface

public ::   sulu_env, init_superlu, prep_superlu, fact_superlu, solv_superlu, free_superlu, close_superlu, Destroy_CompCol_Matrix, Destroy_SuperNode_Matrix, &
            SAMEPATTERN, FACTORED, DOFACT

contains
   !=========================================================================================
   !< @note
   ! **Subroutine to set the default LU behaviour**
   !
   ! + sulu%options%Fact               = ```DOFACT```
   ! + sulu%options%Equil              = ```YES```
   ! + sulu%options%ColPerm            = ```COLAMD```
   ! + sulu%options%DiagPivotThresh    = ```1.0```
   ! + sulu%options%Trans              = ```NOTRANS```
   ! + sulu%options%IterRefine         = ```NOREFINE```
   ! + sulu%options%SymmetricMode      = ```NO```
   ! + sulu%options%PivotGrowth        = ```NO```
   ! + sulu%options%ConditionNumber    = ```NO```
   ! + sulu%options%PrintStat          = ```YES```
   !
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine init_superlu(sulu)
   implicit none
   type(SULU_ENV), intent(inout) :: sulu

      call set_default_options(sulu%options)

   return
   endsubroutine init_superlu


   !=========================================================================================
   !> @note **Subroutine to prepare the [[SULU_ENV]] components**
   !-----------------------------------------------------------------------------------------
   subroutine prep_superlu(sulu)
   implicit none
   type(sulu_env), intent(inout) :: sulu
      integer(kind=I4) :: nn, nz, nb

      nn = sulu%n
      nz = sulu%nz
      nb = sulu%nrhs

      if (.not.allocated( sulu%perm_c )) then
         allocate( sulu%perm_c(1:nn) )
         allocate( sulu%perm_r(1:nn) )
         allocate( sulu%etree( 1:nn) )

         allocate( sulu%RR(1:nn) )
         allocate( sulu%CC(1:nn) )

         allocate( sulu%ferr(1:nb) )
         allocate( sulu%berr(1:nb) )

         allocate(   sulu%rpg(1:nb) )
         allocate( sulu%rcond(1:nb) )

         allocate( sulu%x(1:nn) )
      endif

      sulu%x(1:nn) = 0

      call dCreate_CompCol_Matrix(  A      = sulu%SMA,     & ! out SuperMatrix
                                    m      = sulu%n,       & ! in int
                                    n      = sulu%n,       & ! in int
                                    nnz    = sulu%nz,      & ! in int
                                    nzval  = sulu%a_elt,   & ! in double dimension()
                                    rowind = sulu%irow,    & ! in int dimension()
                                    colptr = sulu%jptr,    & ! in int dimension()
                                    stype  = SLU_NC,       & ! in int
                                    dtype  = SLU_D,        & ! in int
                                    mtype  = SLU_GE        & ! in int
                                 )

      call dCreate_Dense_Matrix( BX    = sulu%smb,  & ! out SuperMatrix
                                 m     = sulu%n,    & ! in int
                                 n     = sulu%nrhs, & ! in int
                                 x     = sulu%b,    & ! in double dimension()
                                 ldx   = sulu%n,    & ! in int
                                 stype = SLU_DN,    & ! in int
                                 dtype = SLU_D,     & ! in int
                                 mtype = SLU_GE     & ! in int
                               )

      call dCreate_Dense_Matrix( BX    = sulu%smx,  & ! out SuperMatrix
                                 m     = sulu%n,    & ! in int
                                 n     = sulu%nrhs, & ! in int
                                 x     = sulu%x,    & ! in double dimension()
                                 ldx   = sulu%n,    & ! in int
                                 stype = SLU_DN,    & ! in int
                                 dtype = SLU_D,     & ! in int
                                 mtype = SLU_GE     & ! in int
                               )

   return
   endsubroutine prep_superlu


   !=========================================================================================
   !< @note
   ! **Subroutine to factorize the system**
   !
   ! note the directives:
   !
   ! + sulu%options%Fact = ```DOFACT```
   ! + sulu%SMB%ncol     = ```0```
   !
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine fact_superlu(sulu, verbose)
   implicit none
   type(SULU_ENV),   intent(inout) :: sulu
   logical(kind=I4), intent(in)    :: verbose

      sulu%lwork = 0
      call StatInit(sulu%stat)

      sulu%options%Fact = DOFACT
      sulu%SMB%ncol     = 0

      call dgssvx(   options            = sulu%options,       & ! superlu_options_t   *options
                     A                  = sulu%SMA,           & ! SuperMatrix         *A
                     perm_c             = sulu%perm_c,        & ! int                 *perm_c
                     perm_r             = sulu%perm_r,        & ! int                 *perm_r
                     etree              = sulu%etree ,        & ! int                 *etree
                     equed              = sulu%equed,         & ! char                *equed
                     R                  = sulu%RR,            & ! double              *R
                     C                  = sulu%CC,            & ! double              *C
                     L                  = sulu%SML,           & ! SuperMatrix         *L
                     U                  = sulu%SMU,           & ! SuperMatrix         *U
                     work               = sulu%work,          & ! void                *work
                     lwork              = sulu%lwork,         & ! int                  lwork
                     B                  = sulu%SMB,           & ! SuperMatrix         *B
                     X                  = sulu%SMX,           & ! SuperMatrix         *X
                     recip_pivot_growth = sulu%rpg,           & ! double              *recip_pivot_growth
                     rcond              = sulu%rcond,         & ! double              *rcond
                     ferr               = sulu%ferr,          & ! double              *ferr
                     berr               = sulu%berr,          & ! double              *berr
                     Glu                = sulu%Glu,           & ! GlobalLU_t          *Glu
                     mem_usage          = sulu%mem_usage,     & ! mem_usage_t         *mem_usage
                     stat               = sulu%stat,          & ! SuperLUStat_t       *stat
                     info               = sulu%info           & ! int                 *info
                  )

      if (verbose) call StatPrint(sulu%stat)
      call StatFree(sulu%stat)

   return
   endsubroutine fact_superlu


   !=========================================================================================
   !< @note
   ! **Subroutine to solve the system**
   !
   ! + If no resolution has yet occured, sulu%first=```true```
   !     * sulu%options%Fact = ```FACTORED```
   !     * sulu%SMB%ncol     = sulu%nrhs (usually ```1```)
   ! + otherwise
   !     * sulu%options%Fact = ```SAMEPATTERN```
   !     * sma, smb and smx are recreated but do not forget that we still have:
   !         - mat%matsulu%irow   => mat%irow
   !         - mat%matsulu%jptr   => mat%jptr
   !         - mat%matsulu%a_elt  => mat%a_elt
   !         - mat%matsulu%b      => mat%b
   !
   ! @endnote
   !
   ! @note
   ! The solution is retrieved with the pointer *store* of type [[NCFORMAT]] which
   ! gives access to [[NCFORMAT:nzval]]
   ! @endnote
   !
   ! @warning
   ! At the end, the memory is released with the dstruction of *sml* and *smu*
   ! @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine solv_superlu(sol_x, sulu, verbose)
   implicit none
   real(kind=R8),    dimension(:), intent(inout) :: sol_x
   type(SULU_ENV),                 intent(inout) :: sulu
   logical(kind=I4),               intent(in)    :: verbose
      type(NCFORMAT), pointer :: Xstore
      real(kind=R8),  pointer :: tabX(:)
      integer(kind=I4)        :: i

      call StatInit(sulu%stat)

      if ( sulu%first ) then
         sulu%options%Fact = FACTORED
         sulu%SMB%ncol     = sulu%nrhs
      else
         sulu%options%Fact = SAMEPATTERN
         call prep_superlu(sulu)
      endif

      call dgssvx(   options            = sulu%options,   & ! superlu_options_t   *options
                     A                  = sulu%SMA,       & ! SuperMatrix         *A
                     perm_c             = sulu%perm_c,    & ! int                 *perm_c
                     perm_r             = sulu%perm_r,    & ! int                 *perm_r
                     etree              = sulu%etree ,    & ! int                 *etree
                     equed              = sulu%equed,     & ! char                *equed
                     R                  = sulu%RR,        & ! double              *R
                     C                  = sulu%CC,        & ! double              *C
                     L                  = sulu%SML,       & ! SuperMatrix         *L
                     U                  = sulu%SMU,       & ! SuperMatrix         *U
                     work               = sulu%work,      & ! void                *work
                     lwork              = sulu%lwork,     & ! int                  lwork
                     B                  = sulu%SMB,       & ! SuperMatrix         *B
                     X                  = sulu%SMX,       & ! SuperMatrix         *X
                     recip_pivot_growth = sulu%rpg,       & ! double              *recip_pivot_growth
                     rcond              = sulu%rcond,     & ! double              *rcond
                     ferr               = sulu%ferr,      & ! double              *ferr
                     berr               = sulu%berr,      & ! double              *berr
                     Glu                = sulu%Glu,       & ! GlobalLU_t          *Glu
                     mem_usage          = sulu%mem_usage, & ! mem_usage_t         *mem_usage
                     stat               = sulu%stat,      & ! SuperLUStat_t       *stat
                     info               = sulu%info       & ! int                 *info
                  )

      call c_f_pointer(sulu%SMX%Store, XStore)
      call c_f_pointer(XStore%nzval, tabX, [XStore%nnz])
      do i = 1, sulu%n
         sol_x(i) = tabX(i)
      enddo
      nullify(Xstore, tabX)

      if (verbose) call StatPrint(sulu%stat)
      call StatFree(sulu%stat)

      call Destroy_SuperNode_Matrix(sulu%SML)
      call Destroy_CompCol_Matrix(  sulu%SMU)

   return
   endsubroutine solv_superlu


   !=========================================================================================
   !< @note Subroutine that actually does nothing yet. Maybe, there will be extra memory that
   ! could be released here?
   !-----------------------------------------------------------------------------------------
   subroutine free_superlu()
   implicit none
   !type(SULU_ENV), intent(inout) :: sulu
   return
   endsubroutine free_superlu


   !=========================================================================================
   !> @note **Subroutine to close the SuperLU process, with memory release**
   !-----------------------------------------------------------------------------------------
   subroutine close_superlu(sulu)
   implicit none
   type(SULU_ENV), intent(inout) :: sulu

      call Destroy_CompCol_Matrix(sulu%SMA)
      call Destroy_Dense_Matrix(  sulu%smb)
      call Destroy_Dense_Matrix(  sulu%smx)

      deallocate( sulu%perm_c )
      deallocate( sulu%perm_r )
      deallocate( sulu%etree  )

      deallocate( sulu%RR )
      deallocate( sulu%CC )

      deallocate( sulu%ferr )
      deallocate( sulu%berr )

      deallocate( sulu%rpg   )
      deallocate( sulu%rcond )

   return
   endsubroutine close_superlu

endmodule sulu_wrapper
