!DEC$ DECLARE

module globvars

!This module contains declarations of global variables
!and constants.

include 'mpif.h'

!Define dprec -- the kind parameter for double precision. This is 
!done in order to simplify possible migration to other platforms 
!where double precision is not the standard 8-byte data type

integer,parameter     :: dprec=8

!=============================================================
!Numerical constants
!=============================================================
real(dprec),parameter :: &
ZERO=0.E0_dprec,     &
ONE=1.E0_dprec,      &
TWO=2.E0_dprec,      &
THREE=3.E0_dprec,    &
FOUR=4.E0_dprec,     &
FIVE=5.E0_dprec,     &
SIX=6.E0_dprec,      &
SEVEN=7.E0_dprec,    &
EIGHT=8.E0_dprec,    &
NINE=9.E0_dprec,     &
TEN=10.0_dprec,      &
ONEHALF=ONE/TWO,     &
ONETHIRD=ONE/THREE,  &
ONEFOURTH=ONE/FOUR,  &
THREEHALF=THREE/TWO, &     
PI=3.14159265358979323846264338327950E0_dprec,     &  
SQRTPI=1.77245385090551602729816748334115E0_dprec, &
Glob_FineStructConst=7.2973525376E-03_dprec          !CODATA 2006 value    

!=============================================================
!Global parameters
!=============================================================

integer,parameter :: Glob_MaxAllowedNumOfPseudoParticles=5
!This is the maximal allowed number of pseudoparticles. If needed,
!this number can be increased. But it is not recommended to use
!a number that is larger than necessary as it may affect the
!performance of the program.

!=============================================================
!Global variables
!These should be set when the program starts 
!=============================================================

integer       Glob_n          !Number of pseudoparticles

!np=n(n+1)/2 - number of independent parameters in a 
!symmetric matrix of size (n x n)   
integer       Glob_np    

!Total number of nonlinear parameters per basis function 
!In case of complex L=0 Gaussians Glob_npt=2*Glob_np
integer       Glob_npt   

!Glob_np_MaxAllowed and Glob_npt_MaxAllowed determine the 
!maximal allowed values for Glob_np and Glob_npt 
integer,parameter :: Glob_np_MaxAllowed=Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2
integer,parameter :: Glob_npt_MaxAllowed=2*Glob_np_MaxAllowed 

real(dprec)   Glob_2raised3n2  !2^(3n/2)

!Glob_CurrBasisSize is a variable whose value equals the current 
!size of the basis
integer       Glob_CurrBasisSize

!Glob_CurrBasisSize is a variable whose value equals the current 
!energy
real(dprec)       Glob_CurrEnergy

!Glob_WhichEigenvalue is a variable that specifies which eigenvalue
!needs to be found. It is used only when the general symmetric
!eigenvalue problem solution method is 'G' 
integer       Glob_WhichEigenvalue

!Glob_ApproxEnergy is a good appriximation to the exact eigenvalue.
!It is used only when the general symmetric eigenvalue problem 
!solution method is 'I'
real(dprec)       Glob_ApproxEnergy

!Glob_EigvalTol is a variable that specifies the accuracy of solving
!the eigenvalue problem. It is used only when the general symmetric
!eigenvalue problem solution method is 'I'
real(dprec)     Glob_EigvalTol

!Glob_InvItParameter is a variable that specifies the factor by which
!an approximate eigenvalue is multiplied when the general symmetric
!eigenvalue problem solution method is 'I'
real(dprec)     Glob_InvItParameter

!Glob_MassMatrix is the mass matrix, M 
real(dprec),allocatable,dimension(:,:),save ::  Glob_MassMatrix

!Glob_Mass is the masses of particles (not pseudoparticles!), M_i
real(dprec),allocatable,dimension(:),save ::  Glob_Mass

!Glob_Mass0 is the total mass of the system, M0
real(dprec)   Glob_Mass0

!Glob_PseudoCharge is the charges of pseudoparticles, qi
real(dprec),allocatable,dimension(:),save ::  Glob_PseudoCharge

!Glob_PseudoCharge0 is the charge of the reference particle, q0
real(dprec)   Glob_PseudoCharge0

!Glob_YOperatorStringLength definess the length of string
!Glob_YOperatorString
integer, parameter  ::  Glob_YOperatorStringLength=100

!Glob_YOperatorString is a string that contains the symbolic 
!expression or the Young operator that is read from the 
!input/output file
character(Glob_YOperatorStringLength)  Glob_YOperatorString

!Glob_LastEigvalTol is a variable that stores the accuracy
!reached last time the eigenvalue problem solver was called.
!It is used only when the general symmetric eigenvalue problem 
!solution method is 'I'
real(dprec)   Glob_LastEigvalTol

!Glob_BestEigvalTol is a variable that stores the best accuracy
!of solving the eigenvalue problem for a certain period of time.
!It is used only when the general symmetric eigenvalue problem 
!solution method is 'I'
real(dprec)   Glob_BestEigvalTol

!Glob_WorstEigvalTol is a variable that stores the best accuracy
!of solving the eigenvalue problem for a certain period of time.
!It is used only when the general symmetric eigenvalue problem 
!solution method is 'I'
real(dprec)   Glob_WorstEigvalTol

!Array Glob_NonlinParam contains the nonlinear parameters of basis 
!functions,
real(dprec),allocatable,dimension(:,:),save :: Glob_NonlinParam

!Array Glob_FuncNum contains the basis function numbers 
integer,allocatable,dimension(:),save :: Glob_FuncNum

!Definition of a new type, Glob_HistoryStep. The data of this
!is used to store the energy and other information obtained
!when the basis size was equal to certain number of functions.

type Glob_HistoryStep
  real(dprec)   Energy
  integer       CyclesDone
  integer       InitFuncAtLastStep
  integer       NumOfEnergyEvalDuringFullOpt
endtype Glob_HistoryStep

!Glob_History is an array of type Glob_HistoryStep that contains
!the information about the energy and basis building and optimization
type(Glob_HistoryStep),allocatable,dimension(:),save :: Glob_History

!4D array Glob_Transposit contains all pair permutation matrices 
!(transpositions). The structure is as follows:
!  Glob_Transposit(1:Glob_n,1:Glob_n,1,2) -- matrix that corresponds to P12
!  Glob_Transposit(1:Glob_n,1:Glob_n,5,5) -- matrix that corresponds to P55 
integer,allocatable,dimension(:,:,:,:),save  :: Glob_Transposit

!Variables Glob_NumYTerms and Glob_NumYHYTerms are the number of 
!independent terms in the Y and Y^{\dagger}Y operator respectively
integer Glob_NumYTerms
integer Glob_NumYHYTerms

!3D arrays Glob_YMatr and Glob_YHYMatr contains all matrices for Y and 
!Y^{\dagger}Y operators. 
!The structure is as follows:
!   Glob_YHYMatr(1:Glob_n,1:Glob_n,5) is the matrix corresponding to the 
!   5-th term of Y^{\dagger}Y operator
integer,allocatable,dimension(:,:,:),save  :: Glob_YMatr
integer,allocatable,dimension(:,:,:),save  :: Glob_YHYMatr

!Arrays Glob_YCoeff and Glob_YHYCoeff contain all coefficients 
!(coefficients of permutations) in the Y and Y^{\dagger}Y operators
integer,allocatable,dimension(:),save  :: Glob_YCoeff
integer,allocatable,dimension(:),save  :: Glob_YHYCoeff

!Variable Glob_GSEPSolutionMethod defines which method
!is currently used to solve the general symmetric eigenvalue
!problem. Initially wet set it to 'U', which means undefined 
character(1)  :: Glob_GSEPSolutionMethod='U'

!Variable Glob_HSLeadDim is the size (leading dimension)
!of H and S matrices, as well as the leading dimension of
!the array that contains derivatives
integer   Glob_HSLeadDim

!Array Glob_H is used to store the Hamiltonian matrix
complex(dprec),allocatable,dimension(:,:),save  :: Glob_H

!Array Glob_S is used to store the overlap matrix
complex(dprec),allocatable,dimension(:,:),save  :: Glob_S

!Array Glob_diagH is used to store the diagonal elements
!of the Hamiltonian matrix 
real(dprec),allocatable,dimension(:),save :: Glob_diagH

!Array Glob_diagS is used to store the diagonal elements
!of the overlap matrix 
real(dprec),allocatable,dimension(:),save :: Glob_diagS

!Array Glob_D is used to store the derivatives of the Hamiltonian
!and the overlap matrix elements. The first index represents
!nonlinear parameters. It ranges from 1 to 4*Glob_np
!if np=Glob_np then
! Glob_D(1:np,i,j)         is dHijdvechLi
! Glob_D(np+1:2*np,i,j)    is dHijdvechBi
! Glob_D(2*np+1:3*np,i,j)  is dSijdvechLi
! Glob_D(3*np+1:4*np,i,j)  is dSijdvechBi 
!If only partial optimization is carried out, that is if only
!last several functions are optimized (= Glob_nfo) then index 2, 
!i, ranges from 1 to Glob_nfo although in fact it should have 
!been equal i+Glob_nfru.  This is done to avoid allocating
!the space which is not used. So
!Glob_D(1:np,i,j) is dH_{i+Glob_nfru,j}/dvechL_{i+Glob_nfru}
! ....
!Index j ranges from 1 to Glob_nfa
complex(dprec),allocatable,dimension(:,:,:),save  :: Glob_D

!Array Glob_c is used to store the eigenvector
complex(dprec),allocatable,dimension(:),save :: Glob_c

!Array Glob_invD is used only when Glob_GSEPSolutionMethod='I'.
!It contains inversed elements of matrix D when
!H-Glob_ApproxEnergy*S is factorized in the Cholessky form, L*D*LH. 
complex(dprec),allocatable,dimension(:),save :: Glob_invD

!Array Glob_LastEigvector contains the eigenvector obtained previous
!time when GSEP was solved using the inverse itaration method. It can
!be used as a initial guess for solving GSEP again, in which case 
!fewer iterations may be needed.
complex(dprec),allocatable,dimension(:),save :: Glob_LastEigvector

!Variable Glob_HSBuffLen defines the length of arrays
!Glob_HBuffer1,Glob_HBuffer2,Glob_SBuffer1, and Glob_SBuffer2
integer  Glob_HSBuffLen

!These 4 arrays are used as buffers to store and send/receive
!the Hamiltonian and the overlap matrix elements
complex(dprec),allocatable,dimension(:),save ::  &
    Glob_HklBuff1,Glob_HklBuff2,Glob_SklBuff1,Glob_SklBuff2

!This array is used as buffer to store and send/receive
!the derivatives of the Hamiltonian and the overlap matrix 
!elements
complex(dprec),allocatable,dimension(:,:),save :: &
    Glob_DkBuff1,Glob_DkBuff2,Glob_DlBuff1,Glob_DlBuff2

!Basis size currently being attempted
integer       Glob_nfa

!Number of functions being added or optimized simultaneously     
integer       Glob_nfo       

!This variable is introduced for the sake of convenience. It should
!be set to be Glob_nfa-Glob_nfo every time any of the 
!latter two change their values. This is basically the number of
!functions to remain unchanged
integer       Glob_nfru  

!Vector Glob_bvc is used for computing particle densities. Its components
!depend on the masses of particles  
real(dprec),allocatable,dimension(:,:),save :: Glob_bvc


!This variable is used to tune routine ZHEGVX (from LAPACK) accuracy
real(dprec)  Glob_AbsTolForZHEGVX

!Variable Glob_MaxIterForGHEPIIS set the limit of iterations allowed
!in GHEPIIS
integer     ::  Glob_MaxIterForGHEPIIS=20

!This parameters defines what is the allowed fraction of
!Eigenvalue problem solution failures for random selection
!process. For example, when Glob_MaxFracOfTrialFailsAllowed=0.15
!then it means no more than 15% of failures allowed. If this
!number is exceeded at some step then the calculations are stopped.
real(dprec),parameter  ::  Glob_MaxFracOfTrialFailsAllowed=0.15*ONE

!This parameter defines how many failures in energy or gradient
!evaluation are allowed during optimization of nonlinear parameters
 integer,parameter   :: Glob_MaxEnergyFailsAllowed=5
! integer,parameter   :: Glob_MaxEnergyFailsAllowed=20

!This parameter defines for how many times basis enlargement 
!routine or the cyclic optimization routine are allowed to repeat 
!random trial and optimization process if the generated 
!function/functions end up being linearly dependent with other 
!functions in the basis (overlap is close to 1.0) or any linear 
!parameters by magnitude are greater then threshold. 
!If this limit is exceeded then the routine stops and shows an 
!error message.
integer,parameter   :: Glob_BadOverlapOrLinCoeffLim=10

!This parameter sets the maximum 2-norm allowed for D (Scaling vector) 
!dot product the very first step that DRMNG attempts. Tuning this value
!is important as it may affect the optimization process very significantly
real(dprec),parameter :: Glob_MaxScStepAllowedInOpt=0.001

!Glob_OverlapPenaltyAllowed is set to be .true. when overlaps are subject
!of penalty for exceeding some threshold
logical           ::  Glob_OverlapPenaltyAllowed=.false.

!Glob_OverlapPenaltyThreshold2 defines the absolute square of the overlap threshold. 
!In case when  Glob_OverlapPenaltyAllowed=.false. and Glob_OverlapPenaltyThreshold2<ONE
!a penalty may be added to the total energy if any of the pair overlaps 
!exceeds Glob_OverlapPenaltyThreshold by magnitude. This global variable
!is changed depending on the user input.
real(dprec)       ::  Glob_OverlapPenaltyThreshold2=ONE

!Glob_MaxOverlapPenalty defines the scale of pair overlap penalties. This global 
!variable is changed depending on the user input. In general, using a value that is
!of the same order of magnitude as the total energy is recommended.
real(dprec)       ::  Glob_MaxOverlapPenalty=ONE

!Glob_TotalOverlapPenalty is equal to the sum of all pair overlap penalties
real(dprec)       ::  Glob_TotalOverlapPenalty=ZERO

!This parameter sets a relative threshold for scaling
!of nonlinear parameters
real(dprec),parameter :: Glob_OptScalingThreshold=1.0E+06_dprec

!This logical parameter defines if the nonlinear parameters of functions
!being optimized are printed on the screen in CycleOptX routines
logical, parameter :: Glob_AreParamPrintedInCycleOptX=.true.

!This parameter defines the level of information printed
!in routines that eliminate basis functions
integer,parameter   :: Glob_ElimRoutPrintSpec=2

!=============================================================
!The following variables are the parameters of the random generator
!that is used to generate new basis function candidates. For the
!description of these see comments in random generator routine,
!GenerateTrialParam
!=============================================================
real(dprec)  Glob_RG_p1
real(dprec)  Glob_RG_s1
real(dprec)  Glob_RG_s2

!=============================================================
!The following variables serve as counters. They are equal the 
!number of times certain routines have been called.
!=============================================================
integer     :: Glob_EnergyGACounter=0
integer     :: Glob_EnergyGBCounter=0
integer     :: Glob_EnergyIACounter=0
integer     :: Glob_EnergyIBCounter=0
!Glob_InvItTempCounter1 and Glob_InvItTempCounter1 are used
!for calculating the average number of inverse iterations
!in certain cycles (when Glob_GSEPSolutionMethod='I')
integer     :: Glob_InvItTempCounter1=0
integer     :: Glob_InvItTempCounter2=0

!=============================================================
!These variables are used to measure time. 
!=============================================================
real(4)                  Glob_TmArray(2)
real(dprec)           :: Glob_TimeSinceStart=ZERO

!=============================================================
!The following arrays and variables provide temporary work space
!for different subroutines
!=============================================================
integer,allocatable,dimension(:)          :: Glob_IntWorkArrForSaveResults
integer                                      Glob_LWorkForZHEGVX
complex(dprec),allocatable,dimension(:)   :: Glob_WorkForZHEGVX
real(dprec),allocatable,dimension(:)      :: Glob_RWorkForZHEGVX
integer,allocatable,dimension(:)          :: Glob_IWorkForZHEGVX
real(dprec),allocatable,dimension(:)      :: Glob_WkGR
complex(dprec),allocatable,dimension(:)   :: Glob_WorkForGHEPIIS
!=============================================================
!Data used to store basis building and optimization program
!=============================================================
!Variable Glob_FileNameLength defines the maximal
!length of file names used throughout code
integer,parameter :: Glob_FileNameLength=40
!The definition of a new type, Glob_BBOPStep
type Glob_BBOPStep
  character(9)  Action
  character(1)  GSEPSolutionMethod
  integer       A
  integer       B
  integer       C
  integer       D
  integer       E
  integer       F
  integer       G
  real(dprec)   Q
  real(dprec)   R
  character(Glob_FileNameLength) FileName1
  character(Glob_FileNameLength) FileName2
  character(Glob_FileNameLength) FileName3
  character(Glob_FileNameLength) FileName4      
endtype Glob_BBOPStep

!This is just a constant string that is used to fill the parameter 
!list in some routines
character(Glob_FileNameLength),parameter :: Glob_FileNameNone='none' 

!Number of Basis Building and Optimization Program Steps
integer    Glob_NumOfBBOPSteps

!Array Glob_BBOP contains the steps of Basis Building
!and Optimization Program
type(Glob_BBOPStep),allocatable,dimension(:),save :: Glob_BBOP

!Current step of BBOP
integer  :: Glob_CurrBBOPStep=0

!=============================================================
!Global variables that define working with files
!=============================================================

!Glob_DataFileName is the name of the input/output file
character(Glob_FileNameLength)          :: Glob_DataFileName='inout.txt'
character(Glob_FileNameLength)          :: Glob_SwapFileName='swapfile.dat'
character(Glob_FileNameLength)          :: Glob_ReallocFileName='realloc.dat'

!These are the variables that specify whether swapping
!is allowed and whether a temporary file should be used
!for array reallocation (the latter reduces memory fragmentation)
logical :: Glob_UseSwapFile=.true.
logical :: Glob_UseReallocFile=.false.

!Parameter Glob_FullOptSaveD controls what data should be
!saved in Hessian file during Full Optimization.
!If Glob_FullOptSaveD=.false. then only the Hessian is saved
!If Glob_FullOptSaveD=.true. then besides the Hessian, the
!program also saves the scaling vector D (used by the optimization 
!routine)
!In the second case, after restart of the full optimization
!the scaling coefficients will be exactly the same as they were
!at the very first start, while when Glob_FullOptSaveD=.false.
!these coefficients are set to some values that depend on the
!current values of nonlinear parameters.
logical,parameter :: Glob_FullOptSaveD=.true.


!=============================================================
!Global variables that contain information about
!identical particles in the system. This information might
!be used to automatically symmetrize the expectation values 
!of two-particle operators. However, an additional symmetrization
!may be required in case when the system under consideration has 
!more complicated symmetry of the Young operator than that 
!dictaded by permutational symmetry of particles 
!=============================================================

!Glob_NumOfIdentPartSets is the number of identical particle sets
!in the system
integer   Glob_NumOfIdentPartSets

!Array Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets) 
!contains the number of particles in each set of identical 
!particles
integer,allocatable,dimension(:),save :: Glob_NumOfPartInIdentPartSet

!Array Glob_IdentPartList contains the list of identical particles
!(their numbers) in each set of identical particles. For example,
!entries Glob_IdentPartList(1:Glob_NumOfPartInIdentPartSet(j),j)
!contain particle numbers that belong to set j.
integer,allocatable,dimension(:,:),save :: Glob_IdentPartList

!Glob_NumOfEqvPairSets is the number of equivalent pseudoparticle 
!pair sets in the system (note that j,j is also considered to be a pair
!even though it involves only pseudoparticle j)
integer   Glob_NumOfNoneqvPairSets

!Array Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets) contains
!the number of pairs in each set of equivalent pairs
integer,allocatable,dimension(:),save :: Glob_NumOfPairsInEqvPairSet

!Array Glob_EqvPairList contains the list of 
!equivalent pairs in each set of equivalent pairs of particles. 
!For example, entries 
!Glob_IdentPartList(1:2,1:Glob_NumOfPairsInEqvPairSet(j),j)
!contain pairs that belong to set j. First index changes from 1 to 2
!designating first and second particle in the pair.
integer,allocatable,dimension(:,:,:),save :: Glob_EqvPairList

!=============================================================
!Global variables used with MPI
!=============================================================

integer     Glob_NumOfProcs !Number of parallel processes
integer     Glob_ProcID     !The ID of a particular process (ranges
!from 0 to Glob_NumOfProcs-1
integer     Glob_MPIErrCode !Error code for MPI routines

!MPI_DPREC -- a new variable for MPI that defines double 
!precision type
integer,parameter    :: MPI_DPREC=MPI_DOUBLE_PRECISION
integer,parameter    :: MPI_DCOMP=MPI_DOUBLE_COMPLEX
integer,parameter    :: MPI_2DPREC=MPI_2DOUBLE_PRECISION

!=============================================================
!Declaring external functions
!=============================================================

real(4), external     :: DTIME
real(dprec), external :: DLAMCH

contains
!=============================================================
!These are routines used exclusively for debugging
!=============================================================
subroutine SaveCMatrixS(n,A,filename)
integer n
complex(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,a1,f8.5,a1,f8.5,a1)',advance='no') '(',real(A(i,j)),',',imag(A(i,j)),')'
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveCMatrixS

subroutine SaveCMatrixSE(n,A,filename)
integer n
complex(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,a2,d12.5,a1,d12.5,a2)',advance='no') ' (',real(A(i,j)),',',imag(A(i,j)),') '
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveCMatrixSE

subroutine SaveCMatrixLE(n,A,filename)
integer n
complex(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,a2,d23.16,a1,d23.16,a2)',advance='no') ' (',real(A(i,j)),',',imag(A(i,j)),') '
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveCMatrixLE

subroutine SaveCVectorS(n,v,filename)
integer n
complex(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,a1,f8.5,a1,f8.5,a1)') '(',real(v(i)),',',imag(v(i)),')'
enddo
close(1)
end subroutine SaveCVectorS

subroutine SaveCVectorSE(n,v,filename)
integer n
complex(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,a2,d12.5,a1,d12.5,a2)') ' (',real(v(i)),',',imag(v(i)),') '
enddo
close(1)
end subroutine SaveCVectorSE

subroutine SaveCVectorLE(n,v,filename)
integer n
complex(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,a2,d23.16,a1,d23.16,a2)') ' (',real(v(i)),',',imag(v(i)),') '
enddo
close(1)
end subroutine SaveCVectorLE

subroutine SaveRMatrixS(n,A,filename)
integer n
real(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,f8.5,1x)',advance='no') A(i,j)
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveRMatrixS

subroutine SaveRMatrixLE(n,A,filename)
integer n
real(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,d23.16,1x)',advance='no') A(i,j)
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveRMatrixLE

subroutine SaveRMatrixSE(n,A,filename)
integer n
real(dprec) A(n,n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  do j=1,n
     write(1,'(1x,d12.5,1x)',advance='no') A(i,j)
  enddo
  write(1,*)
enddo
close(1)
end subroutine SaveRMatrixSE

subroutine SaveRVectorS(n,v,filename)
integer n
real(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,f8.5)') v(i)
enddo
close(1)
end subroutine SaveRVectorS

subroutine SaveRVectorSE(n,v,filename)
integer n
real(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,d12.5)') v(i)
enddo
close(1)
end subroutine SaveRVectorSE

subroutine SaveRVectorLE(n,v,filename)
integer n
real(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,d23.16)') v(i)
enddo
close(1)
end subroutine SaveRVectorLE

subroutine SaveIVector(n,v,filename)
integer n
integer(dprec) v(n)
character(*) filename
integer i,j
open(1,file=filename,status='replace')
do i=1,n
  write(1,'(1x,i6)') v(i)
enddo
close(1)
end subroutine SaveIVector

end module globvars
