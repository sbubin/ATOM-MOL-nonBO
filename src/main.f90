!DEC$ DECLARE
program main

use workproc

!Local variables
integer      i,iw,Kstart,Kstop,Kstep,OpenFileErr,OptimizationType
real(dprec)  r  

!These variables are used to set random number generators
integer RNSeedSize
integer,allocatable :: Seed(:)  

!Initialize MPI
call MPI_INIT(Glob_MPIErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write (*,*) 'Program Expilitly Correlated Complex Gaussians has started'
  write (*,*) 'Number of parallel processes running',Glob_NumOfProcs
  write (*,*)
endif

call ReadIOFile()
call ProgramDataInit()

!Seed the random number generators
call random_seed()
call random_seed(size=RNSeedSize)
allocate(Seed(RNSeedSize))
call random_seed(get=Seed(1:RNSeedSize))
Seed(1:RNSeedSize)=Seed(1:RNSeedSize)+Glob_ProcID
call random_seed(put=Seed(1:RNSeedSize))
call random_number(r)
r=drnor_start(nint(r*20000)+Glob_ProcID)

!Empty swap file as it may contain some garbage left
!after last run (if there was a failure)
if ((Glob_ProcID==0).and.Glob_UseSwapFile) then
  open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
  if (OpenFileErr==0) then
    write(1) 'Swap file is empty'
    close(1)
  endif
endif

OptimizationType=1

do i=1,Glob_NumOfBBOPSteps
  Glob_CurrBBOPStep=i
  Glob_ApproxEnergy=Glob_CurrEnergy*Glob_InvItParameter

  select case (Glob_BBOP(i)%Action)
  
  case('BASIS_ENL')
    Kstart=Glob_BBOP(i)%A
	Kstop=Glob_BBOP(i)%B
    Kstep=Glob_BBOP(i)%C
	if (Kstop>Glob_CurrBasisSize) then
	  if (Kstart<=Glob_CurrBasisSize+1) then
        Kstart=Glob_CurrBasisSize+1
	    select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call BasisEnlG(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
		                 Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R) 
	    case('I')
          call BasisEnlI(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
		                 Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R) 
	    endselect
	  else
        if (Glob_ProcID==0) then
		  write(*,*) 'Error in main: incorrect BBOP step ',i
		  write(*,*) 'One or more parameters in BASIS_ENL are incorrect'
        endif
	  endif
    endif
    
  case('OPT_CYCLE')
    if (Glob_ProcID==0) then
      if ((Glob_BBOP(i)%C>Glob_BBOP(i)%A).or.(Glob_BBOP(i)%B>Glob_BBOP(i)%C)) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
		  write(*,*) 'One or more parameters in OPT_CYCLE are incorrect'
      endif
    endif
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
	    (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
      if (Glob_History(Glob_CurrBasisSize)%CyclesDone<Glob_BBOP(i)%F) then
	    select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')  
	      call OptCycleG(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
		       Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R)    	  
        case('I')
	      call OptCycleI(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
		       Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R) 
	    endselect  	  
      endif
	endif
	
  case('FULL_OPT1')
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
	    (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call FullOpt1G(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
	      real(Glob_BBOP(i)%E,dprec),real(Glob_BBOP(i)%F,dprec),Glob_BBOP(i)%FileName1)	  
      case('I')
	    call FullOpt1I(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
	      real(Glob_BBOP(i)%E,dprec),real(Glob_BBOP(i)%F,dprec),Glob_BBOP(i)%FileName1)	 
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'One or more parameters in FULL_OPT1 are incorrect'
      endif
	endif
	
  case('ELIM_LCFN')
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call EliminateLittleContribFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LCFN'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'One or more parameters in ELIM_LCFN are incorrect'
      endif
	endif

  case('ELIM_LND1')
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call EliminateLinDepFunc1(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LND1'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'One or more parameters in ELIM_LND1 are incorrect'
      endif
	endif

  case('SEPR_LND1')
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call SeparateLinDepFunc1(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'One or more parameters in SEPR_LND1 are incorrect'
      endif
	endif

  case('SEPR_FLCF')
    if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call SeparateFuncLargeCoeff(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'One or more parameters in SEPR_FLCF are incorrect'
      endif
	endif

  case('EXPC_VALS')
    if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call ExpectationValuesG(1,Glob_FileNameNone,Glob_FileNameNone,Glob_FileNameNone,Glob_FileNameNone)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I has not been implemented in EXPC_VALS'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'Second parameters in EXPC_VALS is incorrect'
      endif
	endif
	
  case('DENSITIES')
    if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
	  select case (Glob_BBOP(i)%GSEPSolutionMethod)
      case('G')      
	    call ExpectationValuesG(1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2,Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4)	  
      case('I')
        if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I has not been implemented in DENSITIES'
	  endselect
	else
      if (Glob_ProcID==0) then
        write(*,*) 'Error in main: incorrect BBOP step ',i
		write(*,*) 'Second parameters in DENSITIES is incorrect'
      endif
	endif	
	
  case('SAVE_FILE')
    if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
	  if (Glob_ProcID==0) then
	    iw=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
	    write(*,*)
	    write(*,'(1x,a21,a<iw>,a3)',advance='no') 'Saving basis in file ',Glob_BBOP(i)%FileName1(1:iw),'...'
	    call SaveResults(Filename=Glob_BBOP(i)%FileName1,Sort='no')
	    write(*,'(a5)') ' done'
	    write(*,*)
	  endif  
    endif	
    
  endselect  
enddo

if (Glob_ProcID==0) then
  if (Glob_UseSwapFile) then
    !Empty swap file as it may contain the data saved after
    !the last basis building and optimization program step
    open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
    if (OpenFileErr==0) then
      write(1) 'Swap file is empty'
      close(1)
	  write(*,*)
      write(*,*) 'Swap file has been cleaned up'
    endif
  endif
  write(*,*) ' '
  write(*,*) 'Basis Building and Optimization Program is completed'
  write(*,*) 'Program stopped'
endif

call MPI_FINALIZE(Glob_MPIErrCode)

end program main