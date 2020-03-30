!DEC$ DECLARE

module workproc

use misc
use matform
use matelem
use linalg

!This module contains basic work subroutines

external ILAENV
integer  ILAENV

contains


subroutine ReadIOFile()
!Subroutine ReadIOFile reads data (nonlinear variational
!parameters and other information) from the input/output 
!file whose name is specified by global variable 
!Glob_DataFileName. If there is no such a file in the 
!current directory then the program stops.

!Local variables:
integer        OpenFileErr
real(dprec)    ReadRealA,ReadRealB
real(dprec),allocatable,dimension(:) :: ReadRealArr
integer        ReadInt
integer        WorkInt(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
real(dprec),allocatable,dimension(:) :: WorkBuffReal
integer,allocatable,dimension(:)     :: WorkBuffInt
integer        i,j,Line,j1,j2,j3,j4
character(70)  ReadChar
logical        ErrorInDataFile,IsBBOPStep

ErrorInDataFile=.false.

if (Glob_ProcID==0) then
  open(1,file=Glob_DataFileName,status='old',iostat=OpenFileErr)
  if (OpenFileErr/=0) then
    print *,'Error in DataFileInit: data file not found - ',Glob_DataFileName
    ErrorInDataFile=.true.
  endif
endif

if (ErrorInDataFile) stop

!Reading information
if (Glob_ProcID==0) Line=0
if (Glob_ProcID==0) then
  write(*,*) 'Reading initial conditions from data file ',Glob_DataFileName
  read(1,*) ReadChar(1:9),ReadInt
  write(*,'(1x,a9,i3)') ReadChar(1:9),ReadInt
  Line=Line+1
  Glob_n=ReadInt-1 !Glob_n is the number of pseudoparticles
  if ((Glob_n<1).or.(ReadChar(1:9)/='PARTICLES')) then
    print *,'Error in data file, line ',Line   
    ErrorInDataFile=.true.
  endif
endif 
if (Glob_n>Glob_MaxAllowedNumOfPseudoParticles) then
  if (Glob_ProcID==0) then
    write (*,*) 'The version of the code you are running was compiled for the case'
    write (*,*) 'when the number of pseudoparticles is smaller or equal to', Glob_MaxAllowedNumOfPseudoParticles
    write (*,*) 'while in the number of pseudoparticles in the system specified by'
    write (*,*) 'the input file is',Glob_n
    write (*,*) 'Please make appropriate changes. Program will now stop.'
  endif
  ErrorInDataFile=.true.
endif
if (ErrorInDataFile) stop
call MPI_BCAST(Glob_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
Glob_np=Glob_n*(Glob_n+1)/2
Glob_npt=Glob_np+Glob_np
Glob_2raised3n2=TWO**((3*Glob_n)/TWO)

allocate(Glob_Mass(Glob_n+1))
if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:6),Glob_Mass(1:Glob_n+1)
  write(*,*) ReadChar(1:6)
  do i=1,Glob_n+1
    write(*,'(1x,d23.16)') Glob_Mass(i)
  enddo
  Line=Line+1
endif
call MPI_BCAST(Glob_Mass,Glob_n+1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

allocate(Glob_PseudoCharge(Glob_n))
if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:7),Glob_PseudoCharge0,Glob_PseudoCharge(1:Glob_n)
  write(*,*) ReadChar(1:7)
  write(*,'(1x,d23.16)') Glob_PseudoCharge0
  do i=1,Glob_n
    write(*,'(1x,d23.16)') Glob_PseudoCharge(i)
  enddo
  Line=Line+1
endif
call MPI_BCAST(Glob_PseudoCharge0,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_PseudoCharge,Glob_n,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:8),Glob_YOperatorString
  i=len_trim(Glob_YOperatorString)
  write(*,'(1x,a8,1x,a<i>)')  & 
             ReadChar(1:8),trim(Glob_YOperatorString)
  Line=Line+1
endif
do i=1,Glob_YOperatorStringLength
   WorkInt(i)=ichar(Glob_YOperatorString(i:i))
enddo
call MPI_BCAST(WorkInt,Glob_YOperatorStringLength,MPI_INTEGER,0, &
               MPI_COMM_WORLD,Glob_MPIErrCode)
do i=1,Glob_YOperatorStringLength
   Glob_YOperatorString(i:i)=char(WorkInt(i))
enddo

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:10),Glob_CurrBasisSize
  write(*,'(1x,a10,1x,i6)')  ReadChar(1:10),Glob_CurrBasisSize
  Line=Line+1
endif
call MPI_BCAST(Glob_CurrBasisSize,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:14),Glob_CurrEnergy
  write(*,'(1x,a14,3x,d23.16)')  ReadChar(1:14),Glob_CurrEnergy
  Line=Line+1
endif
call MPI_BCAST(Glob_CurrEnergy,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:16),Glob_WhichEigenvalue
  write(*,'(1x,a16,1x,i6)')  ReadChar(1:16),Glob_WhichEigenvalue
  Line=Line+1
endif
call MPI_BCAST(Glob_WhichEigenvalue,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:16),Glob_EigvalTol
  write(*,'(1x,a16,1x,d23.16)')  ReadChar(1:16),Glob_EigvalTol
  Line=Line+1
endif
call MPI_BCAST(Glob_EigvalTol,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:14),Glob_InvItParameter
  write(*,'(1x,a14,3x,d23.16)')  ReadChar(1:14),Glob_InvItParameter
  Line=Line+1
endif
call MPI_BCAST(Glob_InvItParameter,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:15),Glob_LastEigvalTol
  write(*,'(1x,a15,2x,d23.16)')  ReadChar(1:15),Glob_LastEigvalTol
  Line=Line+1
endif
call MPI_BCAST(Glob_LastEigvalTol,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:15),Glob_BestEigvalTol
  write(*,'(1x,a15,2x,d23.16)')  ReadChar(1:15),Glob_BestEigvalTol
  Line=Line+1
endif
call MPI_BCAST(Glob_BestEigvalTol,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:16),Glob_WorstEigvalTol
  write(*,'(1x,a16,1x,d23.16)')  ReadChar(1:16),Glob_WorstEigvalTol
  Line=Line+1
endif
call MPI_BCAST(Glob_WorstEigvalTol,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:15),Glob_RG_p1,Glob_RG_s1,Glob_RG_s2
  write(*,'(1x,a15)') ReadChar(1:15)
  write(*,'(1x,d23.16)') Glob_RG_p1
  write(*,'(1x,d23.16)') Glob_RG_s1
  write(*,'(1x,d23.16)') Glob_RG_s2
  Line=Line+1
endif
call MPI_BCAST(Glob_RG_p1,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RG_s1,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RG_s2,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,'(a70)')   ReadChar(1:70)
  write(*,'(1x,a70)')  ReadChar(1:70)
  read(1,'(a70)')   ReadChar(1:70)
  write(*,'(1x,a70)')  ReadChar(1:70)
  read(1,'(a70)')   ReadChar(1:70)
  write(*,'(1x,a70)')  ReadChar(1:70)
  Line=Line+3
endif

!Reading Basis Building and Optimization Program
ReadChar(1:70)=' '
if (Glob_ProcID==0) then
  Glob_NumOfBBOPSteps=0
  IsBBOPStep=.true.
  do while (IsBBOPStep)
    read(1,*) ReadChar(1:9)
	if ((ReadChar(1:9)=='BASIS_ENL').or.(ReadChar(1:9)=='OPT_CYCLE').or.  &
	    (ReadChar(1:9)=='FULL_OPT1').or.(ReadChar(1:9)=='EXPC_VALS').or.  &
		(ReadChar(1:9)=='ELIM_LCFN').or.(ReadChar(1:9)=='ELIM_LND1').or.  &
		(ReadChar(1:9)=='SEPR_LND1').or.(ReadChar(1:9)=='SEPR_FLCF').or.  &
		(ReadChar(1:9)=='DENSITIES').or.(ReadChar(1:9)=='SAVE_FILE')) then
	  Glob_NumOfBBOPSteps=Glob_NumOfBBOPSteps+1
    else
      IsBBOPStep=.false.
	endif
  enddo
  do i=1,Glob_NumOfBBOPSteps+1 
    backspace 1
  enddo
endif
call MPI_BCAST(Glob_NumOfBBOPSteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
allocate(Glob_BBOP(Glob_NumOfBBOPSteps))

if (Glob_ProcID==0) then
  do i=1,Glob_NumOfBBOPSteps
    read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod    
  enddo
  do i=1,Glob_NumOfBBOPSteps
    backspace 1
  enddo
  do i=1,Glob_NumOfBBOPSteps
    select case (Glob_BBOP(i)%Action(1:9))
    case('BASIS_ENL')
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod , &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C, &
			Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
      write(*,'(1x,a9,1x,a1,5(1x,i6),1x,d23.16,1x,d23.16)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
    case('OPT_CYCLE')
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C, &
			Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G, &
			Glob_BBOP(i)%Q,Glob_BBOP(i)%R
      write(*,'(1x,a9,1x,a1,7(1x,i6),1x,d23.16,1x,d23.16)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E, &
			Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
    case('FULL_OPT1')			
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D, &
            Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
			Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,a1,4(1x,i6),2(1x,d23.16),2(1x,i6),1x,a<j>)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
			Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%FileName1(1:j)			
    case('EXPC_VALS')
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
      write(*,'(1x,a9,1x,a1,1(1x,i6))') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A	  
    case('DENSITIES')
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName3(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)
      j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
      j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
      j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))                              	
      write(*,'(1x,a9,1x,a1,1x,i6,1x,a<j1>,1x,a<j2>,1x,a<j3>,1x,a<j4>)')              &
            Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,  &
	        Glob_BBOP(i)%FileName1(1:j1),Glob_BBOP(i)%FileName2(1:j2),                &	        
	        Glob_BBOP(i)%FileName3(1:j3),Glob_BBOP(i)%FileName4(1:j4)	              
	case('ELIM_LCFN') 
	  read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &		
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,a1,1x,i6,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%FileName1(1:j)
	case('ELIM_LND1') 
	  read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &		
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,a1,1x,i6,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%FileName1(1:j)
	case('SEPR_LND1') 
	  read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &		
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,a1,1x,i6,d23.16,1x,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1(1:j)
    case('SEPR_FLCF')
	  read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &		
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,a1,1x,i6,d23.16,1x,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1(1:j)
    case('SAVE_FILE')
      read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)   
      j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(*,'(1x,a9,1x,i6,1x,a<j>)') Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:j)                 
	endselect			     
  enddo
  read(1,'(a70)')   ReadChar(1:70)
  write(*,'(1x,a70)')  ReadChar(1:70)
endif
do i=1,Glob_NumOfBBOPSteps
  do j=1,9
    WorkInt(j)=ichar(Glob_BBOP(i)%Action(j:j))
  enddo
  call MPI_BCAST(WorkInt,9,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,9
     Glob_BBOP(i)%Action(j:j)=char(WorkInt(j))
  enddo
  j=ichar(Glob_BBOP(i)%GSEPSolutionMethod)
  call MPI_BCAST(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  Glob_BBOP(i)%GSEPSolutionMethod=char(j)
  call MPI_BCAST(Glob_BBOP(i)%A,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%B,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%C,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%D,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%E,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%F,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%G,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_BBOP(i)%Q,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  call MPI_BCAST(Glob_BBOP(i)%R,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_BBOP(i)%FileName1(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_BBOP(i)%FileName1(j:j)=char(WorkInt(j))
  enddo
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_BBOP(i)%FileName2(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_BBOP(i)%FileName2(j:j)=char(WorkInt(j))
  enddo  
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_BBOP(i)%FileName3(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_BBOP(i)%FileName3(j:j)=char(WorkInt(j))
  enddo  
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_BBOP(i)%FileName4(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_BBOP(i)%FileName4(j:j)=char(WorkInt(j))
  enddo  
enddo

allocate(Glob_History(Glob_CurrBasisSize))
allocate(Glob_FuncNum(Glob_CurrBasisSize))
allocate(Glob_NonlinParam(Glob_npt,Glob_CurrBasisSize))

if (Glob_CurrBasisSize==0) then
    !Stop reading if basis size is zero
	close(1)
	return
endif

allocate(WorkBuffReal(Glob_CurrBasisSize))
allocate(WorkBuffInt(Glob_CurrBasisSize))

!Reading history
if (Glob_ProcID==0) then
  do i=1,Glob_CurrBasisSize
    read(1,*)  j,Glob_History(i)%Energy,Glob_History(i)%CyclesDone, &
	      Glob_History(i)%InitFuncAtLastStep, &
		  Glob_History(i)%NumOfEnergyEvalDuringFullOpt
  enddo
  read(1,'(a70)')   ReadChar(1:70)
  write(*,'(1x)')
endif
do i=1,Glob_CurrBasisSize
  WorkBuffReal(i)=Glob_History(i)%Energy
enddo
call MPI_BCAST(WorkBuffReal,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
do i=1,Glob_CurrBasisSize
  Glob_History(i)%Energy=WorkBuffReal(i)
enddo
do i=1,Glob_CurrBasisSize
  WorkBuffInt(i)=Glob_History(i)%CyclesDone
enddo
call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
do i=1,Glob_CurrBasisSize
  Glob_History(i)%CyclesDone=WorkBuffInt(i)
enddo
do i=1,Glob_CurrBasisSize
  WorkBuffInt(i)=Glob_History(i)%InitFuncAtLastStep
enddo
call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
do i=1,Glob_CurrBasisSize
  Glob_History(i)%InitFuncAtLastStep=WorkBuffInt(i)
enddo
do i=1,Glob_CurrBasisSize
  WorkBuffInt(i)=Glob_History(i)%NumOfEnergyEvalDuringFullOpt
enddo
call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
do i=1,Glob_CurrBasisSize
  Glob_History(i)%NumOfEnergyEvalDuringFullOpt=WorkBuffInt(i)
enddo

deallocate(WorkBuffReal)
deallocate(WorkBuffInt)

!Reading the basis functions
if (Glob_ProcID==0) then
  do i=1,Glob_CurrBasisSize
    read(1,*) j,Glob_NonlinParam(1:Glob_npt,i)
  enddo
endif
call MPI_BCAST(Glob_NonlinParam,Glob_npt*Glob_CurrBasisSize, &
                MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 

!Setting the function numbers as they were read
do i=1,Glob_CurrBasisSize
  Glob_FuncNum(i)=i
enddo

close(1)

end subroutine ReadIOFile


subroutine SaveResults(Filename,Sort)
!Subroutine SaveResults saves current results of the
!calculations in the input/output file.
!  Input parameters :
!  FileName - Optional file name, where the results need to be saved. If it is
!         not passed then the results are saved in the file whose name is defined
!         by global variable Glob_DataFileName.
!  Sort - Optional string variable. If its value is equal to 'yes' or 'YES', 
!         then the basis functions are sorted out in according to 
!         their numbers. In this case a work array is needed. The 
!         subroutine uses array Glob_IntWorkArrForSaveResults. Thus,
!         in order to do calls with Sort='yes' one first needs to
!         allocate that global array. The length of this integer 
!         array should be at least Glob_CurrBasisSize. If the value of
!         Sort is not equal to 'yes' or 'YES', or it is simply not passed, then
!         the basis function will not be sorted out. In this case the allocation
!         of global array Glob_IntWorkArrForSaveResults is not necessary.
!Arguments:
character(*)  FileName,Sort
optional   :: FileName,Sort
!Local variables: 
integer i,j,j1,j2,j3,j4
logical SortNeeded

if (Glob_ProcID==0) then
  if (present(FileName)) then
    open(1,file=FileName,status='replace')
  else
    open(1,file=Glob_DataFileName,status='replace')
  endif
  write(1,'(1x,a9,i3)') 'PARTICLES',Glob_n+1
  write(1,'(1x,a6,<Glob_n+1>(1x,d23.16))') 'MASSES',Glob_Mass(1:Glob_n+1)
  write(1,'(1x,a7,1x,d23.16,<Glob_n>(1x,d23.16))') 'CHARGES',  &
              Glob_PseudoCharge0,Glob_PseudoCharge(1:Glob_n)
  i=len_trim(Glob_YOperatorString)
  write(1,'(1x,a8,1x,a<i>)')    'SYMMETRY',Glob_YOperatorString(1:i)
  write(1,'(1x,a10,1x,i6)')     'BASIS_SIZE',Glob_CurrBasisSize
  write(1,'(1x,a14,3x,d23.16)') 'CURRENT_ENERGY',Glob_CurrEnergy
  write(1,'(1x,a16,1x,i6)')     'WHICH_EIGENVALUE',Glob_WhichEigenvalue
  write(1,'(1x,a16,1x,d23.16)') 'EIGVAL_TOLERANCE',Glob_EigvalTol
  write(1,'(1x,a14,3x,d23.16)') 'INVITPARAMETER',Glob_InvItParameter
  write(1,'(1x,a15,2x,d23.16)') 'LAST_EIGVAL_TOL',Glob_LastEigvalTol
  write(1,'(1x,a15,2x,d23.16)') 'BEST_EIGVAL_TOL',Glob_BestEigvalTol
  write(1,'(1x,a16,1x,d23.16)') 'WORST_EIGVAL_TOL',Glob_WorstEigvalTol
  write(1,'(1x,a15,1x,d23.16,1x,d23.16,1x,d23.16)') &
      'GENERATOR_PARAM',Glob_RG_p1,Glob_RG_s1,Glob_RG_s2
  write(1,'(1x,a30)') '=============================='
  i=Glob_CurrBasisSize
  write(1,'(1x,i6,1x,d23.16,1x,i6,1x,i6,1x,i7)')  i,  &
       Glob_History(i)%Energy,Glob_History(i)%CyclesDone,  &
	   Glob_History(i)%InitFuncAtLastStep,  &
	   Glob_History(i)%NumOfEnergyEvalDuringFullOpt
  write(1,'(1x,a30)') '=============================='
  do i=1,Glob_NumOfBBOPSteps
    select case (Glob_BBOP(i)%Action(1:9))
    case('BASIS_ENL')
      write(1,'(1x,a9,1x,a1,5(1x,i6),1x,d23.16,1x,d23.16)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
    case('OPT_CYCLE')
      write(1,'(1x,a9,1x,a1,7(1x,i6),1x,d23.16,1x,d23.16)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E, &
			Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
    case('FULL_OPT1')
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,a1,4(1x,i6),2(1x,d23.16),2(1x,i6),1x,a<j>)') Glob_BBOP(i)%Action(1:9), &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
			Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
			Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%FileName1(1:j) 
    case('EXPC_VALS')
      write(1,'(1x,a9,1x,a1,1x,i6)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A	
    case('DENSITIES')
	  j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))  
	  j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength)) 
	  j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength)) 
	  j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)) 	  	  	    
      write(1,'(1x,a9,1x,a1,1x,i6,1x,a<j1>,1x,a<j2>,1x,a<j3>,1x,a<j4>)')  &
            Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,     &
            Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:j1),                  &
            Glob_BBOP(i)%FileName2(1:j2),Glob_BBOP(i)%FileName3(1:j3),    &	        
	        Glob_BBOP(i)%FileName4(1:j4) 
	case('ELIM_LCFN') 
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,a1,1x,i6,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%FileName1(1:j)
	case('ELIM_LND1') 
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,a1,1x,i6,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%FileName1(1:j)
	case('SEPR_LND1') 
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,a1,1x,i6,d23.16,1x,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1(1:j)
    case('SEPR_FLCF')
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,a1,1x,i6,d23.16,1x,d23.16,1x,a<j>)') Glob_BBOP(i)%Action(1:9),  &
	        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%Q, &
            Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1(1:j)
    case('SAVE_FILE')
	  j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
      write(1,'(1x,a9,1x,i6,1x,a<j>)') Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:j)            
	endselect			     
  enddo
  write(1,'(1x,a30)') '=============================='
  do i=1,Glob_CurrBasisSize
    write(1,'(1x,i6,1x,d23.16,1x,i6,1x,i6,1x,i7)') i,Glob_History(i)%Energy, &
	    Glob_History(i)%CyclesDone,Glob_History(i)%InitFuncAtLastStep, &
		Glob_History(i)%NumOfEnergyEvalDuringFullOpt
  enddo
  write(1,'(1x,a30)') '=============================='
  SortNeeded=.false.
  if (present(Sort)) then
    if ((Sort=='yes').or.(Sort=='YES')) SortNeeded=.true.
  endif
  if (SortNeeded) then 
    do i=1,Glob_CurrBasisSize
      Glob_IntWorkArrForSaveResults(Glob_FuncNum(i))=i
    enddo
    do i=1,Glob_CurrBasisSize
      write(1,'(1x,i6,<Glob_npt>(1x,d23.16))') i,  &
	     Glob_NonlinParam(1:Glob_npt,Glob_IntWorkArrForSaveResults(i))
    enddo
  else
    do i=1,Glob_CurrBasisSize
      write(1,'(1x,i6,<Glob_npt>(1x,d23.16))') i,Glob_NonlinParam(1:Glob_npt,i)
    enddo
  endif
  close(1)
endif

end subroutine SaveResults



subroutine ProgramDataInit()
!Subroutine ProgramDataInit initializes some data needed for
!calculations. It should be called in the beginning of the
!calculations, right after reading input/output file.

integer       n,npart
integer       i,j,k,p,q,t,s,w,ii,jj,kk
character(1)  c1,cc1
real(dprec)   Mtot
integer       StrLen,NumFactY
integer       TotNumOfYTerms,TotNumOfYHYTerms,CurrNumOfTerms
integer       L,R,FirstLPos,LastRPos,MaxNumTermsInFact
integer,allocatable,dimension(:)      :: TempSymCoeff,TempSymCoeff1,NumTermsInYOpFact
integer,allocatable,dimension(:,:,:)  :: TempSymMatr,TempSymMatr1
integer,allocatable,dimension(:,:)    :: Matr1,Matr2,Matr3,Matr4
integer                               :: Coeff,Cf3
character(Glob_YOperatorStringLength),allocatable,dimension(:) :: YOpStr, YHOpStr
integer       pi,pj,pt,ps
logical       AreTermsIdentical
integer,allocatable,dimension(:)      :: IdentParticleSet
integer,allocatable,dimension(:,:)    :: IdentPseudoPartPairSet

if (Glob_ProcID==0) write(*,*) 'Initializing program data'

!Initialization of some machine-dependent numbers
Glob_AbsTolForZHEGVX=2*DLAMCH('S')

n=Glob_n
npart=n+1

!Constructing Glob_MassMatrix
allocate(Glob_MassMatrix(n,n))
Glob_MassMatrix(1:n,1:n)=ONEHALF/Glob_Mass(1)
do i=1,n
  Glob_MassMatrix(i,i)=Glob_MassMatrix(i,i)+ONEHALF/Glob_Mass(i+1)
enddo

!Determine the components of vector Glob_bvc
!that is used in evaluation of particle densities
allocate(Glob_bvc(n,npart)) 
Mtot=sum(Glob_Mass(1:npart))
do i=1,npart
  Glob_bvc(1:n,i)=-Glob_Mass(2:n+1)/Mtot
enddo
do i=2,npart
  Glob_bvc(i-1,i)=Glob_bvc(i-1,i)+ONE
enddo

!write(*,*)
!write(*,*) 'Mtot=',Mtot
!write(*,*)
!do i=1,npart
!  do j=1,n
!    write(*,'(1x,a9,i1,a1,i1,a2,d23.16)') 'Glob_bvc(',j,',',i,')=',Glob_bvc(j,i)
!  enddo
!  write(*,*)
!enddo
!write(*,*)

!Constructing all Pij transposition matrices 
!
!  P1i  (i/=1) has the following view:
!
!                    i-1
!
!      |  1   0   0  -1   0   0   0  |
!      |                             |
!      |  0   1   0  -1   0   0   0  |
!      |                             |
!      |  0   0   1  -1   0   0   0  |            
!      |                             |          This particular
!P1i = |  0   0   0  -1   0   0   0  |  i-1     matrix is P15 for the 
!      |                             |          case of 8 particles
!      |  0   0   0  -1   1   0   0  |
!      |                             |
!      |  0   0   0  -1   0   1   0  |
!      |                             |
!      |  0   0   0  -1   0   0   1  |
!
!
!  Pij has the following form:
!
!                i-1         j-1
!
!      |  1   0   0   0   0   0   0  |
!      |                             |
!      |  0   1   0   0   0   0   0  |
!      |                             |
!      |  0   0   0   0   0   1   0  |  i-1
!      |                             |
!Pij = |  0   0   0   1   0   0   0  |          This particular  
!      |                             |          matrix is P47 for the
!      |  0   0   0   0   1   0   0  |          case of 8 particles
!      |                             |
!      |  0   0   1   0   0   0   0  |  j-1
!      |                             |
!      |  0   0   0   0   0   0   1  |
!
allocate(Glob_Transposit(n,n,npart,npart))
!First set all of them to be unit matrices
Glob_Transposit(1:n,1:n,1:npart,1:npart)=0
do i=1,npart
  do j=1,npart
    do k=1,n
      Glob_Transposit(k,k,i,j)=1	  
	enddo
  enddo
enddo
!Now continue depending on type of transposition (P1i or Pij)
do i=2,npart
  Glob_Transposit(1:n,i-1,1,i)=-1
enddo
do i=2,npart
  do j=i+1,npart
    Glob_Transposit(i-1,i-1,i,j)=0
	Glob_Transposit(j-1,j-1,i,j)=0
	Glob_Transposit(j-1,i-1,i,j)=1
	Glob_Transposit(i-1,j-1,i,j)=1
    Glob_Transposit(i-1,i-1,j,i)=0
	Glob_Transposit(j-1,j-1,j,i)=0
	Glob_Transposit(j-1,i-1,j,i)=1
	Glob_Transposit(i-1,j-1,j,i)=1
  enddo
enddo

!Constructing the Young operator based on the content of
!a string variable Glob_YOperatorString

!First we throw away spaces and multiplication signs
!from Glob_YOperatorString
StrLen=len_trim(Glob_YOperatorString)
do i=1,StrLen
  c1=Glob_YOperatorString(i:i)
  if ((c1==' ').or.(c1=='*')) then
    do j=i,StrLen-1
      Glob_YOperatorString(j:j)=Glob_YOperatorString(j+1:j+1)
	enddo
    Glob_YOperatorString(j:j)=' '
  endif
enddo
StrLen=len_trim(Glob_YOperatorString)

!Checking for wrong symbols in Glob_YOperatorString
do i=1,StrLen
  c1=Glob_YOperatorString(i:i)
  if ((c1/='1').and.(c1/='2').and.(c1/='3').and.(c1/='4').and.(c1/='5').and. &
      (c1/='6').and.(c1/='7').and.(c1/='8').and.(c1/='9').and.(c1/='0').and. &
	  (c1/='P').and.(c1/='+').and.(c1/='-').and.(c1/='*').and.(c1/=')').and. &
	  (c1/='(')) then
    write(*,*) 'Error in ProgramDataInit: the Young operator expression'
    write(*,*) 'contains wrong symbols'
    stop
  endif
enddo

!Checking if the number of left and right brackets is the same
!and counting how many brackets there are
L=0
R=0
do i=1,StrLen
  if (Glob_YOperatorString(i:i)==')') R=R+1
  if (Glob_YOperatorString(i:i)=='(') L=L+1
enddo
if (R/=L) then
  write(*,*) 'Error in ProgramDataInit: the numer of left and right brackets in the'
  write(*,*) 'Young operator is different'
  stop
endif

!NumFactY is the number of factors in the Young operator,
!FirstLPos is the position of the first left bracket,
!LastRPos is the position of the last right bracket,
if (R/=0) then
  FirstLPos=scan(Glob_YOperatorString(1:StrLen),'(')
  LastRPos=scan(Glob_YOperatorString(1:StrLen),')',back=.true.)
  NumFactY=R
  i=0
  do j=1,R
	k=0
	i=i+1
    c1=Glob_YOperatorString(i:i)
    do while (c1/='(')
	  if (c1/='*') k=1
	  i=i+1
	  c1=Glob_YOperatorString(i:i)
	enddo 
	if (k==1) NumFactY=NumFactY+1
	do while (c1/=')')
	  i=i+1
	  c1=Glob_YOperatorString(i:i)
	enddo
  enddo
  if (Glob_YOperatorString(StrLen:StrLen)/=')') NumFactY=NumFactY+1
else
  NumFactY=1 
endif

!Splitting Glob_YOperatorString into an array of smaller
!strings, YOpStr. Each column of this array will contain just
!one factor, with no brackets. A '+' or a '-' sign is added in 
!front of the first term in a factor if needed. Multiplication
!signs are dropped .
allocate(YOpStr(NumFactY))
do i=1,NumFactY
  YOpStr(i)=' '
enddo 
if (R==0) then
  c1=Glob_YOperatorString(1:1)
  if ((c1/='+').or.(c1/='-')) then
    YOpStr(1)(1:1)='+'
	YOpStr(1)(2:StrLen+1)=Glob_YOperatorString(1:StrLen)
  else
	YOpStr(1)(1:StrLen)=Glob_YOperatorString(1:StrLen) 
  endif
else
  i=1
  k=1
  p=i
  q=0
  c1=Glob_YOperatorString(i:i)
  if ((c1/='(').and.(c1/='+').and.(c1/='-')) then
    q=1
    YOpStr(k)(1:1)='+'
  endif
  do while (Glob_YOperatorString(i:i)/='(')
    i=i+1
  enddo
  if (i>1) then
    YOpStr(k)(p+q:i-1+q)=Glob_YOperatorString(p:i-1)
    !if (YOpStr(k)(i-1+q:i-1+q)=='*') YOpStr(k)(i-1+q:i-1+q)=' '
    k=k+1
  endif
  do j=1,R
    i=i+1
	p=i
	q=0
    c1=Glob_YOperatorString(i:i)
    if ((c1/=')').and.(c1/='+').and.(c1/='-')) then
      q=1
      YOpStr(k)(1:1)='+'
	endif
	do while (Glob_YOperatorString(i:i)/=')')
      i=i+1
	enddo
    YOpStr(k)(1+q:i+q-p)=Glob_YOperatorString(p:i-1)
    k=k+1
	i=i+1
    c1=Glob_YOperatorString(i:i)
	if (c1=='*') then
      i=i+1
      c1=Glob_YOperatorString(i:i)
	endif
	if ((c1/='(').and.(c1/=' '))  then
	  p=i
      YOpStr(k)(1:1)='+'	
	  do while ((c1/='(').and.(c1/=' '))
	    i=i+1
        c1=Glob_YOperatorString(i:i)
      enddo
      YOpStr(k)(2:i+1-p)=Glob_YOperatorString(p:i-1)
	  k=k+1
    endif
  enddo
endif

!Print all factors in the Young operator
!j=StrLen+1
!do i=1,NumFactY
!  write (*,'(1x,i3,1x,a3,a<j>)') i,':  ',YOpStr(i)(1:j)
!enddo

!Creating an array that contains all the factors of the 
!Y^{\dagger} operator. Basically, Y^{\dagger} is the reversed Y (i.e.
!the order of all factors is reversed as well as 
!permutation products (if there are any) in each factor come
!in reverse order.
allocate(YHOpStr(NumFactY))
do i=1,NumFactY
  YHOpStr(i)=' '
enddo 
do i=NumFactY,1,-1
  s=NumFactY-i+1
  j=1
  c1=YOpStr(s)(j:j)
  do while (c1/=' ')
    if (c1=='P') then
	  k=0 !k counts the number of Permutations in the current term
	  t=0
      do while ((c1/='+').and.(c1/='-').and.(c1/=' '))
        if (c1=='P') k=k+1
		t=t+1
		c1=YOpStr(s)(j+t:j+t)
	  enddo
	  do t=1,k
        YHOpStr(i)(j+3*(k-t):j+3*(k-t)+2)=YOpStr(s)(j+3*(t-1):j+3*(t-1)+2)
	  enddo
	  j=j+3*k
	else
      YHOpStr(i)(j:j)=c1
      j=j+1
	endif
	c1=YOpStr(s)(j:j)
  enddo
enddo

!Counting how many terms there are in each factor of the Young
!operator, as well as the total number of terms in the nonsimplified
!Young operator
allocate(NumTermsInYOpFact(NumFactY))
TotNumOfYTerms=1
do k=1,NumFactY
  j=0
  do i=1,Glob_YOperatorStringLength
    if ((YOpStr(k)(i:i)=='+').or.(YOpStr(k)(i:i)=='-')) j=j+1
  enddo
  NumTermsInYOpFact(k)=j
  TotNumOfYTerms=TotNumOfYTerms*j
  TotNumOfYHYTerms=TotNumOfYTerms*TotNumOfYTerms
enddo
if (Glob_ProcID==0) then
  write(*,'(1x,a60,i7)')  'Total number of terms in the nonsimplified Y operator:      ',TotNumOfYTerms
  write(*,'(1x,a60,i7)')  'Total number of terms in the nonsimplified Y^{+}Y operator: ',TotNumOfYHYTerms
endif

allocate(Matr1(1:n,1:n))
allocate(Matr2(1:n,1:n))
allocate(Matr3(1:n,1:n))
allocate(Matr4(1:n,1:n))

!Multiplying all factors in YOpStr and placing actual matrices
!and coefficients in arrays Glob_YMatr and Glob_YCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.  
CurrNumOfTerms=NumTermsInYOpFact(NumFactY)
allocate(TempSymCoeff(CurrNumOfTerms)) 
allocate(TempSymMatr(n,n,CurrNumOfTerms))  
do j=NumFactY,1,-1
  !reading the current factor
  k=0
  i=1
  c1=YOpStr(j)(i:i)
  p=i
  do while (c1/=' ')
    i=i+1
    c1=YOpStr(j)(i:i)
    do while ((c1/='P').and.(c1/='+').and.(c1/='-'))
      i=i+1
      c1=YOpStr(j)(i:i)
    enddo
    if (i-p>1) then
      read(YOpStr(j)(p:i-1),*) Coeff
    else
      if (YOpStr(j)(i-1:i-1)=='+') then
        Coeff=1
	  else
        Coeff=-1
	  endif
    endif 
    Matr1=Glob_Transposit(1:n,1:n,1,1)
    do while (c1=='P')
      read(YOpStr(j)(i+1:i+1),*) p
      read(YOpStr(j)(i+2:i+2),*) q
	  Matr2(1:n,1:n)=Glob_Transposit(1:n,1:n,p,q)
	  Matr4(1:n,1:n)=Matr1(1:n,1:n)
	  do ii=1,n
	    do jj=1,n
	      w=0
	      do kk=1,n
	        w=w+Matr2(ii,kk)*Matr4(kk,jj)
	      enddo
	      Matr1(ii,jj)=w
	    enddo
	  enddo
	  i=i+3
      c1=YOpStr(j)(i:i)
    enddo
    k=k+1
    p=i
    if (j/=NumFactY) then
      if (k==1) then
	    Matr3(1:n,1:n)=Matr1(1:n,1:n)
        Cf3=Coeff
	  else
	    do s=1,t
          Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
          q=t*(k-1)+s
	      do ii=1,n
	        do jj=1,n
	          w=0
	          do kk=1,n
	            w=w+Matr2(ii,kk)*Matr1(kk,jj)
	          enddo
	          TempSymMatr(ii,jj,q)=w
	        enddo
	      enddo
	      TempSymCoeff(q)=Coeff*TempSymCoeff(s)
	    enddo
      endif
	else
      TempSymMatr(1:n,1:n,k)=Matr1(1:n,1:n)
      TempSymCoeff(k)=Coeff
    endif
  enddo
  if (j/=NumFactY) then
    do s=1,t
      Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
	  do ii=1,n
	     do jj=1,n
	        w=0
	        do kk=1,n
	          w=w+Matr2(ii,kk)*Matr3(kk,jj)
	        enddo
	        TempSymMatr(ii,jj,s)=w
	     enddo
	  enddo      
	  TempSymCoeff(s)=Cf3*TempSymCoeff(s)
    enddo
  endif
  !mark the identical terms (adding their coefficients
  !and setting all of them but one to zero)
  t=CurrNumOfTerms
  do i=1,CurrNumOfTerms
    if (TempSymCoeff(i)==0) cycle
    do s=i+1,CurrNumOfTerms
      if (TempSymCoeff(s)==0) cycle    
      if (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
        TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
        if (TempSymCoeff(i)==0) t=t-1
	    TempSymCoeff(s)=0
	    t=t-1
	  endif
    enddo
  enddo     
  !reallocate arrays containing symmetry terms
  !to allow for multiplication by the next factor  
  if (j/=1) then   
    allocate(TempSymCoeff1(t)) 
    allocate(TempSymMatr1(n,n,t)) 
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        TempSymCoeff1(s)=TempSymCoeff(i)
        TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      endif
    enddo
    CurrNumOfTerms=t*NumTermsInYOpFact(j-1)
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    allocate(TempSymCoeff(CurrNumOfTerms)) 
    allocate(TempSymMatr(n,n,CurrNumOfTerms))     
    TempSymCoeff(1:t)=TempSymCoeff1(1:t)
    TempSymMatr(1:n,1:n,1:t)=TempSymMatr1(1:n,1:n,1:t)
    deallocate(TempSymCoeff1)
    deallocate(TempSymMatr1)     
  endif    
enddo

Glob_NumYTerms=t
allocate(Glob_YCoeff(Glob_NumYTerms)) 
allocate(Glob_YMatr(n,n,Glob_NumYTerms))
s=0
do i=1,CurrNumOfTerms
  if (TempSymCoeff(i)/=0) then
    s=s+1
    Glob_YCoeff(s)=TempSymCoeff(i)
    Glob_YMatr(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
  endif
enddo
deallocate(TempSymCoeff)
deallocate(TempSymMatr)
if (Glob_ProcID==0) then
  write(*,'(1x,a60,i7)')  'Total number of terms in the simplified Y operator:         ',Glob_NumYTerms
endif

!Now doing the same thing for Y^{\dagger}Y operator, that
!is expanding it and collecting identical terms

!Multiplying all factors in YHOpStr by already existing
!matrices and coefficients of Y. and placing actual matrices
!and coefficients in arrays Glob_YHYMatr and Glob_YHYCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.  
CurrNumOfTerms=NumTermsInYOpFact(1)*Glob_NumYTerms
allocate(TempSymCoeff(CurrNumOfTerms)) 
allocate(TempSymMatr(n,n,CurrNumOfTerms))  
!TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
!TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)
TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)
t=Glob_NumYTerms
do j=NumFactY,1,-1
  !reading the current factor
  k=0
  i=1
  c1=YHOpStr(j)(i:i)
  p=i
  do while (c1/=' ')
    i=i+1
    c1=YHOpStr(j)(i:i)
    do while ((c1/='P').and.(c1/='+').and.(c1/='-'))
      i=i+1
      c1=YHOpStr(j)(i:i)
    enddo
    if (i-p>1) then
      read(YHOpStr(j)(p:i-1),*) Coeff
    else
      if (YHOpStr(j)(i-1:i-1)=='+') then
        Coeff=1
	  else
        Coeff=-1
	  endif
    endif 
    Matr1=Glob_Transposit(1:n,1:n,1,1)
    do while (c1=='P')
      read(YHOpStr(j)(i+1:i+1),*) p
      read(YHOpStr(j)(i+2:i+2),*) q
	  Matr2(1:n,1:n)=Glob_Transposit(1:n,1:n,p,q)
	  Matr4(1:n,1:n)=Matr1(1:n,1:n)
	  do ii=1,n
	    do jj=1,n
	      w=0
	      do kk=1,n
	        w=w+Matr2(ii,kk)*Matr4(kk,jj)
	      enddo
	      Matr1(ii,jj)=w
	    enddo
	  enddo
	  i=i+3
      c1=YHOpStr(j)(i:i)
    enddo
    k=k+1
    p=i
    if (k==1) then
	  Matr3(1:n,1:n)=Matr1(1:n,1:n)
      Cf3=Coeff
	else
	  do s=1,t
        Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
        q=t*(k-1)+s
	    do ii=1,n
	      do jj=1,n
	          w=0
	          do kk=1,n
	            w=w+Matr2(ii,kk)*Matr1(kk,jj)
	          enddo
	          TempSymMatr(ii,jj,q)=w
	      enddo
	    enddo
	    TempSymCoeff(q)=Coeff*TempSymCoeff(s)
	  enddo
    endif
  enddo
  do s=1,t
    Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
	do ii=1,n
	   do jj=1,n
	      w=0
	      do kk=1,n
	        w=w+Matr2(ii,kk)*Matr3(kk,jj)
	      enddo
	      TempSymMatr(ii,jj,s)=w
	   enddo
	enddo        
	TempSymCoeff(s)=Cf3*TempSymCoeff(s)
  enddo
  !mark the identical terms (adding their coefficients
  !and setting all of them but one to zero)
  t=CurrNumOfTerms
  do i=1,CurrNumOfTerms
    if (TempSymCoeff(i)==0) cycle
    do s=i+1,CurrNumOfTerms
      if (TempSymCoeff(s)==0) cycle    
      if (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
        TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
        if (TempSymCoeff(i)==0) t=t-1
	    TempSymCoeff(s)=0
	    t=t-1
	  endif
    enddo
  enddo     
  !reallocate arrays containing symmetry terms
  !to allow for multiplication by the next factor
  if (j/=1) then   
    allocate(TempSymCoeff1(t)) 
    allocate(TempSymMatr1(n,n,t)) 
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        TempSymCoeff1(s)=TempSymCoeff(i)
        TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      endif
    enddo
    CurrNumOfTerms=t*NumTermsInYOpFact(NumFactY-j+2)
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    allocate(TempSymCoeff(CurrNumOfTerms)) 
    allocate(TempSymMatr(n,n,CurrNumOfTerms))     
    TempSymCoeff(1:t)=TempSymCoeff1(1:t)
    TempSymMatr(1:n,1:n,1:t)=TempSymMatr1(1:n,1:n,1:t)
    deallocate(TempSymCoeff1)
    deallocate(TempSymMatr1)     
  endif  
enddo

Glob_NumYHYTerms=t
allocate(Glob_YHYCoeff(Glob_NumYHYTerms)) 
allocate(Glob_YHYMatr(n,n,Glob_NumYHYTerms))
s=0
do i=1,CurrNumOfTerms
  if (TempSymCoeff(i)/=0) then
    s=s+1
    Glob_YHYCoeff(s)=TempSymCoeff(i)
    Glob_YHYMatr(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
  endif
enddo
deallocate(TempSymCoeff)
deallocate(TempSymMatr)
if (Glob_ProcID==0) then
  write(*,'(1x,a60,i7)')  'Total number of terms in the simplified Y^{+}Y operator:    ',Glob_NumYHYTerms
endif

!Print all independent Y operator matrices and coefficients
!open(1,file='symterms_new.txt',status='replace')
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYTerms
!    write(1,'(10x,a2,i3,a4,i3)') 'i=',i,'  C=',Glob_YCoeff(i)
!    do j=1,n
!      write(1,'(<n>(1x,i2))') Glob_YMatr(j,1:n,i)
!    enddo 
!    write(1,*)
!  enddo
!endif
!write(1,*)
!write(1,*) '-----'
!Print all independent Y^{\dagger}Y operator matrices and coefficients
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYHYTerms
!    write(1,'(10x,a2,i3,a4,i3)') 'i=',i,'  C=',Glob_YHYCoeff(i)
!    do j=1,n
!      write(1,'(<n>(1x,i2))') Glob_YHYMatr(j,1:n,i)
!    enddo 
!    write(1,*)
!  enddo
!endif
!close(1)
!stop

!open(1,file='symterms_new.txt',status='replace')
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYHYTerms
!    write(1,'(1x,a14,i2,a2,i3)') 'Glob_YHYCoeff(',i,')=',Glob_YHYCoeff(i)
!  enddo
!  do i=1,Glob_NumYHYTerms
!    do j=1,n
!      do k=1,n
!        write(1,'(1x,a13,i2,a1,i2,a1,i2,a2,i2)') 'Glob_YHYMatr(', k , ',' , j , ',' , i, ')=', Glob_YHYMatr(k,j,i)
!      enddo
!    enddo
!  enddo
!endif
!close(1)
!stop

deallocate(Matr4)
deallocate(Matr3)
deallocate(Matr2)
deallocate(Matr1)
deallocate(NumTermsInYOpFact)
deallocate(YHOpStr)
deallocate(YOpStr)

!Now we determine which particles are identical. This determination
!is based on the input values of masses and charges only. The information
!about the sets of identical particles may be needed for proper
!symmetrization of expectation values of operators that involve
!two-paricle quantities (such as interparticle distances).
!The set of particles to which particle i belongs is labelled by IdentParticleSet(i)
!If IdentParticleSet(i)=IdentParticleSet(j) then it means that
!particles i and j are identical. The largest value in IdentParticleSet
!gives the total number of identical particle sets 
allocate(IdentParticleSet(npart))
IdentParticleSet(1)=1
k=1
do i=2,npart
  s=0
  j=0
  do while ((j<i-1).and.(s==0))
    j=j+1
    if (j>1) then
      if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge(j-1)==Glob_PseudoCharge(i-1))) then
        IdentParticleSet(i)=IdentParticleSet(j)
		s=1
	  endif
    else 
	  !j=1 case
      if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge0==Glob_PseudoCharge(i-1))) then
        IdentParticleSet(i)=IdentParticleSet(j)
		s=1
	  endif
	endif
  enddo	
  if (s==0) then
    k=k+1
    IdentParticleSet(i)=k
  endif
enddo
Glob_NumOfIdentPartSets=maxval(IdentParticleSet(1:npart))

!Below we determine which pairs of pseudoparticles are identical.
!The information about this is stored in array IdentPseudoPartPairSet(1:n,1:n)
!Diagonal elements do not actually designate pairs of pseudoparticles but
!rather a single pseudoparticle, which corresponds to a certain pair of particles. 
!If IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(k,l) then
!it means these pairs ij and kl should be equivalent.
!The largest value of array IdentPseudoPartPairSet gives the number
!of nonequivalent pairs.
allocate(IdentPseudoPartPairSet(1:n,1:n))
IdentPseudoPartPairSet(1:n,1:n)=0
k=0
do i=1,n
  do j=i,n   
    if (i==j) then
      pi=1; pj=j+1
	else
      pi=i+1; pj=j+1
	endif
    w=0
    do s=1,i
	  if (s==i) then
		q=j-1
      else
        q=n
	  endif
      do t=s,q
	    if (w==1) cycle
        if (s==t) then
          ps=1; pt=t+1
		else
          ps=s+1; pt=t+1
		endif
        if ((IdentParticleSet(ps)==IdentParticleSet(pi)).and. &
		    (IdentParticleSet(pt)==IdentParticleSet(pj))) then
          w=1
          IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(s,t)
        endif
	  enddo
	enddo
	if (w==0) then
      k=k+1
      IdentPseudoPartPairSet(i,j)=k
	endif
	enddo
enddo
Glob_NumOfNoneqvPairSets=maxval(IdentPseudoPartPairSet(1:n,1:n))

!Now we create arrays Glob_NumOfPartInIdentPartSet
!and Glob_IdentPartList
allocate(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
allocate(Glob_IdentPartList(npart,Glob_NumOfIdentPartSets))
Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)=0
Glob_IdentPartList(1:npart,1:Glob_NumOfIdentPartSets)=0
do i=1,npart
  k=IdentParticleSet(i)
  Glob_NumOfPartInIdentPartSet(k)=Glob_NumOfPartInIdentPartSet(k)+1
  Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k),k)=i
enddo

!Create arrays Glob_NumOfPairsInEqvPairSet and
!Glob_EqvPairList
allocate(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
allocate(Glob_EqvPairList(2,n*(n+1)/2,Glob_NumOfNoneqvPairSets))
Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets)=0
Glob_EqvPairList(1:2,1:n*(n+1)/2,1:Glob_NumOfNoneqvPairSets)=0
do i=1,n
  do j=i,n
    k=IdentPseudoPartPairSet(i,j)
    Glob_NumOfPairsInEqvPairSet(k)=Glob_NumOfPairsInEqvPairSet(k)+1
    Glob_EqvPairList(1,Glob_NumOfPairsInEqvPairSet(k),k)=i
    Glob_EqvPairList(2,Glob_NumOfPairsInEqvPairSet(k),k)=j
  enddo
enddo

deallocate(IdentPseudoPartPairSet)
deallocate(IdentParticleSet)

end subroutine ProgramDataInit



subroutine GenerateTrialParam(nfun,x,k,method_used)
!Subroutine GenerateTrialParam generates nonlinear 
!parameters for trial functions. It uses three global variables
!to do that: Glob_RG_p1, Glob_RG_s1, Glob_RG_s2. 
!Glob_RG_p1 is the probability of using the first generation method, 
!so 1-Glob_RG_p1 is the probability of using the second one. Both 
!methods select a random function from the existing basis (or several 
!functions if nfun>1). Then, the first method generates nonlinear 
!parameters that are independently normally distributed around the 
!corresponding parameters of the selected function with standard 
!deviation Glob_RG_s1*x(i), where x(i) is the i-th nonlinear parameter 
!of the selected function. The second method generates nonlinear 
!parameters in the following way: it generates a random number, r, from 
!a normal (0,Glob_RG_s2) distribution and then multiplies all x(i) by (1+r). 
!It also makes sure that 0.8<|1+r|>1.2 as we do not want almost linearly
!dependent functions in the basis.
!  Input parameter:
!   nfun - number of trial functions whose parameters
!          are to be generated; 
!  Output parameters :
!    x(1:Glob_npt,1:nfun) - a 2D-array containing generated
!          nonlinear parameters of the trial functions;
!      k - specifies which function in the existing 
!          basis was selected to generate a candidate. In case when
!          nfun>1 k specifies which was the first function (out of 
!          a set) selected. In case Glob_CurrBasisSize==0 k is set to 0.
!    method_used - specifies which generation method was used (1 or 2).
!          In case Glob_CurrBasisSize==0 method_used is set to 0.

!Arguments :
integer        nfun
real(dprec)    x(Glob_npt,nfun)
integer        k, method_used       

!Local variables :
integer        i,j,p
real(dprec)    r,sumf
!Constants that define the uniform distribution
!in the case of Glob_CurrBasisSize==0 
real(dprec)  :: Lmin=-0.5
real(dprec)  :: Lmax= 0.5
real(dprec)  :: Bmin=-0.5
real(dprec)  :: Bmax= 0.5


if (Glob_CurrBasisSize==0) then !making uniform distribution
  do i=1,nfun
    do j=1,Glob_np
      call random_number(r)
	  x(j,i)=r*(Lmax-Lmin)+Lmin
      call random_number(r)
	  x(Glob_np+j,i)=r*(Bmax-Bmin)+Bmin	  
    enddo
  enddo
  method_used=0
  k=0
else
  if (Glob_CurrBasisSize<nfun) then
    call random_number(r)
	if (r<GLob_RG_p1) then
      !method 1
      do i=1,nfun
	    call random_number(r)
	    k=int(r*(Glob_CurrBasisSize))+1
		do j=1,Glob_npt
           x(j,i)=(Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j,k)
		enddo
      enddo
      method_used=1
    else
      !method 2
	  r=Glob_RG_s2*drnor()+ONE
	  do while ((abs(r)>0.8E0_dprec).and.(abs(r)<1.2))
        r=Glob_RG_s2*drnor()+ONE
	  enddo
      do i=1,nfun
	    call random_number(r)
	    k=int(r*(Glob_CurrBasisSize))+1
		do j=1,Glob_npt
           x(j,i)=r*Glob_NonlinParam(j,k)
		enddo
      enddo
      method_used=2
	endif
  else
	call random_number(r)
	k=int(r*(Glob_CurrBasisSize-nfun+1))+1
    call random_number(r)
    if (r<GLob_RG_p1) then
      !method 1
      do i=1,nfun
		do j=1,Glob_npt
           x(j,i)=(Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j,k+i-1)
		enddo
      enddo
	  method_used=1
	else
      !method 2
	  r=Glob_RG_s2*drnor()+ONE
	  do while ((abs(r)>0.8E0_dprec).and.(abs(r)<1.2))
        r=Glob_RG_s2*drnor()+ONE
	  enddo
      do i=1,nfun
		do j=1,Glob_npt
           x(j,i)=r*Glob_NonlinParam(j,k+i-1)
		enddo
      enddo
      method_used=2
	endif
  endif
endif

end subroutine GenerateTrialParam



subroutine ComputeOverlapPenalty(MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty)
!Subroutine ComputeOverlapPenalty computes overlap penalty. Overlap
!penalty may be used during full optimization of all nonlinear parameters 
!in order to prevent severe pair linear dependencies of basis functions. 
!The penalty function is P = sum_{i=1,K; j=Glob_nfru+1,K; i<j} Pij, where 
!  Pij=(abs(Sij)^2-t^2)*b/(1-t^2),  Sij^2>t^2
!  Pij=0, Sij>=t^2
!Here 
!  b=MaxPairOverlapPenalty
!  t^2=OverlapThreshold2
!  Sij is the pair overlap value
!The result is returned in TotalPenalty
!Arguments:
real(dprec)  MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty
!Local variables
integer      i,j,k,nbands,leftover
real(dprec)  pen_coeff,tp,aux
logical      oddband

tp=ZERO
pen_coeff=MaxPairOverlapPenalty/(ONE-OverlapThreshold2)
if (Glob_NumOfProcs==1) then
  !In case of a single MPI process we do not split work
  do i=1,Glob_nfru
    do j=Glob_nfru+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) tp=tp+pen_coeff*(aux-OverlapThreshold2)
    enddo
  enddo  
  do i=Glob_nfru+1,Glob_nfa
    do j=i+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) tp=tp+pen_coeff*(aux-OverlapThreshold2)    
    enddo
  enddo      
else
  !In case of more than one MPI process we split the work
  !more or less evenly between processes
  do i=1+Glob_ProcID,Glob_nfru,Glob_NumOfProcs
    do j=Glob_nfru+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) tp=tp+pen_coeff*(aux-OverlapThreshold2)    
    enddo
  enddo
  nbands=Glob_nfo/Glob_NumOfProcs
  leftover=mod(Glob_nfo,Glob_NumOfProcs)
  oddband=.true.
  do k=1,nbands
    if (oddband) then
      i=Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
      oddband=.false.
    else
      i=Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
      oddband=.true.
    endif 
    do j=i+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) tp=tp+pen_coeff*(aux-OverlapThreshold2)
    enddo
  enddo !k
  if (leftover>1) then
    if (oddband) then
      i=Glob_nfa-Glob_ProcID+1
    else
      i=Glob_nfa-leftover+1+Glob_ProcID
    endif
    if ((i<Glob_nfa).and.(i>Glob_nfa-leftover)) then
      do j=i+1,Glob_nfa
        aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
        if (aux>OverlapThreshold2) tp=tp+pen_coeff*(aux-OverlapThreshold2)
      enddo
    endif        
  endif
endif    
call MPI_ALLREDUCE(tp,TotalPenalty,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

end subroutine ComputeOverlapPenalty



subroutine ComputeOverlapPenaltyAndAddGradient(MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty,WkGR)
!Subroutine ComputeOverlapPenalty computes overlap penalty as well as the 
!addition to the energy gradient due to the penalty. Overlap
!penalty may be used during full optimization of all nonlinear parameters 
!in order to prevent severe pair linear dependencies of basis functions. 
!The penalty function is P = sum_{i=1,K; j=Glob_nfru+1,K; i<j} Pij, where 
!  Pij=(abs(Sij)^2-t^2)*b/(1-t^2),  Sij^2>t^2
!  Pij=0, Sij>=t^2
!Here 
!  b=MaxPairOverlapPenalty
!  t^2=OverlapThreshold2
!  Sij is the pair overlap value
!The overlap penalty is returned in TotalPenalty. The gradient is added to WkGR.
!Note that even though the subroutine divides the work between MPI processes, it
!does not combine the results for the gradient. It is assumed that such an operation will be done
!manually right after ComputeOverlapPenaltyAndAddGradient is called. This is done to avoid
!combining the gradient components twice (once after the gradient is computed and once
!after the gradient addition due to the overalp penalties is computed). 
!Arguments:
real(dprec)  MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty,WkGR(:)
!Local variables
integer      i,j,k,m,nbands,leftover
real(dprec)  pen_coeff,tp,aux,temp1,temp2
logical      oddband

tp=ZERO
pen_coeff=MaxPairOverlapPenalty/(ONE-OverlapThreshold2)
temp1=TWO*pen_coeff
if (Glob_NumOfProcs==1) then
  !In case of a single MPI process we do not split work
  do i=1,Glob_nfru
    do j=Glob_nfru+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) then
        tp=tp+pen_coeff*(aux-OverlapThreshold2)
        temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
        do m=1,Glob_npt
          WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*conjg(  Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2  )  )
        enddo        
      endif  
    enddo
  enddo  
  do i=Glob_nfru+1,Glob_nfa
    do j=i+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) then
        tp=tp+pen_coeff*(aux-OverlapThreshold2)
        temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
        do m=1,Glob_npt
          WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*(  Glob_D(Glob_npt+m,i-Glob_nfru,j)-ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2  )  )
        enddo   
        do m=1,Glob_npt
          WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*conjg(  Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2  )  )
        enddo        
      endif     
    enddo
  enddo    
else
  !In case of more than one MPI process we split the work
  !more or less evenly between processes
  do i=1+Glob_ProcID,Glob_nfru,Glob_NumOfProcs
    do j=Glob_nfru+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) then
        tp=tp+pen_coeff*(aux-OverlapThreshold2)
        temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
        do m=1,Glob_npt
          WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*conjg(  Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2  )  )
        enddo        
      endif     
    enddo
  enddo
  nbands=Glob_nfo/Glob_NumOfProcs
  leftover=mod(Glob_nfo,Glob_NumOfProcs)
  oddband=.true.
  do k=1,nbands
    if (oddband) then
      i=Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
      oddband=.false.
    else
      i=Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
      oddband=.true.
    endif 
    do j=i+1,Glob_nfa
      aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
      if (aux>OverlapThreshold2) then
        tp=tp+pen_coeff*(aux-OverlapThreshold2)
        temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
        do m=1,Glob_npt
          WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*(  Glob_D(Glob_npt+m,i-Glob_nfru,j)-ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2  )  )
        enddo     
        do m=1,Glob_npt
          WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*conjg(  Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2  )  )        
        enddo   
      endif  
    enddo
  enddo !k
  if (leftover>1) then
    if (oddband) then
      i=Glob_nfa-Glob_ProcID+1
    else
      i=Glob_nfa-leftover+1+Glob_ProcID
    endif
    if ((i<Glob_nfa).and.(i>Glob_nfa-leftover)) then
      do j=i+1,Glob_nfa
        aux=real(Glob_S(j,i))*real(Glob_S(j,i))+imag(Glob_S(j,i))*imag(Glob_S(j,i))
        if (aux>OverlapThreshold2) then
          tp=tp+pen_coeff*(aux-OverlapThreshold2)
          temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
          do m=1,Glob_npt
            WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*(  Glob_D(Glob_npt+m,i-Glob_nfru,j)-ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2  )  )
          enddo     
          do m=1,Glob_npt
            WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m)+temp1*real(  Glob_S(j,i)*conjg(  Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2  )  )        
          enddo        
        endif  
      enddo
    endif        
  endif
endif    
call MPI_ALLREDUCE(tp,TotalPenalty,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

end subroutine ComputeOverlapPenaltyAndAddGradient



function EnergyGA(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyGA computes the energy of the system under 
!consideration. Routine ZHEGVX from LAPACK is used to solve 
!GSEP. The function finds the energy level defined by the 
!global variable Glob_WhichEigenvalue. The calculations are
!done with a basis of Nmax functions. It is assumed that 
!nonlinear parameters of the basis functions are stored in 
!global array Glob_NonlinParam and all the matrix elements of the
!first Nmin-1 functions have already been calculated and stored in 
!proper global arrays. If no matrix elements have been calculated,
!one should set Nmin=0. In case when all matrix elements have been
!computed calculated and it is only necessary to solve GSEP one 
!needs to set AreMatElemNeeded=.false.
!Upon successfull exit ErrorCode should be equal to 0. If
!ErrorCode = Nmax + i then the leading minor of order i of S is 
!not positive definite.

real(dprec)    EnergyGA
!Arguments:
integer        Nmin,Nmax
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables:
integer        i,j
real(dprec)    Evalue, EVs(1)
complex(dprec) Z(1)
integer        NumOfEigvalsFound
integer        IFAIL     



if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
if (Nmax==1) then 
  EnergyGA=Glob_diagH(1)
  ErrorCode=0
else
  !Copying H and S matrix elements from the lower triangles 
  !to the upper ones. The diagonals are copied from the global
  !arrays where they are stored
  do i=1,Nmax
    do j=1,i-1
	  Glob_H(j,i)=conjg(Glob_H(i,j))
	enddo
	Glob_H(i,i)=Glob_diagH(i)
  enddo
  do i=1,Nmax
    do j=1,i-1
	  Glob_S(j,i)=conjg(Glob_S(i,j))
	enddo
	Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
  enddo
  if (Glob_ProcID==0) then
	call ZHEGVX(1,'N','I','U',Nmax,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
       ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	   NumOfEigvalsFound,EVs,Z,Nmax,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
       Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
     !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     !$                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
     !$                   LWORK, RWORK, IWORK, IFAIL, INFO )   
    Evalue=EVs(1)  
  endif
  if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_OverlapPenaltyAllowed) then
    EnergyGA=Evalue+Glob_TotalOverlapPenalty   
  else  
    EnergyGA=Evalue
  endif 
endif

Glob_EnergyGACounter=Glob_EnergyGACounter+1

end function EnergyGA



function EnergyGAM(Nmin,Nmax,AreMatElemNeeded,ErrorCode)  
!Function EnergyGAM is essentially the same as EnergyGA.
!The only difference is that EnergyGAM also computes the
!linear coefficients (stored in Glob_c).

real(dprec)    EnergyGAM
!Arguments:
integer        Nmin,Nmax
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables:
integer        i,j
real(dprec)    Evalue, EVs(1)
integer        NumOfEigvalsFound
integer        IFAIL


if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
if (Nmax==1) then 
  EnergyGAM=Glob_diagH(1)
  Glob_c(1)=ONE
  ErrorCode=0
else
  !Copying H and S matrix elements from the lower triangles 
  !to the upper ones. The diagonals are copied from the global
  !arrays where they are stored
  do i=1,Nmax
    do j=1,i-1
	  Glob_H(j,i)=conjg(Glob_H(i,j))
	enddo
	Glob_H(i,i)=Glob_diagH(i)
  enddo
  do i=1,Nmax
    do j=1,i-1
	  Glob_S(j,i)=conjg(Glob_S(i,j))
	enddo
	Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
  enddo
  if (Glob_ProcID==0) then
	call ZHEGVX(1,'V','I','U',Nmax,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &   
       ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	   NumOfEigvalsFound,EVs,Glob_c,Nmax,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &   
       Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
     !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     !$                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
     !$                   LWORK, RWORK, IWORK, IFAIL, INFO ) 
    Evalue=EVs(1)  
  endif
  if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Glob_c,Nmax,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  if (Glob_OverlapPenaltyAllowed) then
    EnergyGAM=Evalue+Glob_TotalOverlapPenalty   
  else  
    EnergyGAM=Evalue
  endif 
endif

Glob_EnergyGACounter=Glob_EnergyGACounter+1

end function EnergyGAM  



subroutine EnergyGB(Evalue,Gradient,AreMatElemNeeded,ErrorCode)
!subroutine EnergyGB computes the energy value and the
!gradient with respect to the nonlinear parameters of the last 
!Glob_nfo functions of the system under consideration. Routine
!ZHEGVX from LAPACK is used to solve GSEP. The subroutine finds
!the energy level defined by the  global variable 
!Glob_WhichEigenvalue, as well as the corresponding gradient. It
!is assumed that nonlinear parameters of the basis functions are
!stored in global array Glob_NonlinParam and all the matrix
!elements of the first Glob_nfru functions have been computed
!and stored in proper global arrays. If no matrix elements have
!been calculated, one should first call ComputeMatElemAndDeriv
!and calculate the matrix elements of the first Glob_nfru
!functions. In case when all matrix elements have been computed
!and it is only necessary to solve GSEP and find gradient, one
!needs to set AreMatElemNeeded=.false.
!Upon successfull exit ErrorCode should be equal to 0. If
!ErrorCode <= Glob_nfa then the eigenvector failed to
!converge. If ErrorCode = Glob_nfa + i then the 
!leading minor of order i of S is not positive definite.
!The data in Gradient is ordered in the following manner:
!Gradient=(dEdvechLi,dEdvechBi,dEdvechL{i+1},...,dEdvechLj,dEdvechBj)
!where i=Glob_nfru+1, j=Glob_nfa

!Arguments:
real(dprec)    Evalue
real(dprec)    Gradient(Glob_npt*Glob_nfo)
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables
integer        nfo,nfa,nfru,npt
integer        i,j,k,l,m,nbands,leftover
logical        oddband
integer        N,NumOfEigvalsFound,IFAIL(1)
real(dprec)    EVs(1)
complex(dprec) W(Glob_npt_MaxAllowed),t
real(dprec)    pen_coeff

nfo=Glob_nfo
nfa=Glob_nfa
npt=Glob_npt
nfru=Glob_nfru

if (AreMatElemNeeded) call ComputeMatElemAndDeriv(nfru+1,nfa)

if (nfa==1) then
  Evalue=Glob_diagH(1)
  Glob_c(1)=ONE
  ErrorCode=0
else
  !Copyinng H and S matrix elements from the lower triangles 
  !to the upper ones. The diagonals are copied from the global
  !arrays where they are stored
  do i=1,nfa
    do j=1,i-1
	  Glob_H(j,i)=conjg(Glob_H(i,j))
	enddo
	Glob_H(i,i)=Glob_diagH(i)
  enddo
  do i=1,nfa
    do j=1,i-1
	  Glob_S(j,i)=conjg(Glob_S(i,j))
	enddo
	Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
  enddo

  if (Glob_ProcID==0) then
	call ZHEGVX(1,'V','I','U',nfa,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
       ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	   NumOfEigvalsFound,EVs,Glob_c,nfa,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
       Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
     !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     !$                   VL, VU, IL, IU, ABSTOL, 
	 !$                   M, W, Z, LDZ, WORK, LWORK, 
     !$                   RWORK, IWORK, IFAIL, INFO ) 
    Evalue=EVs(1)
  endif
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Glob_c,nfa,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

!Computing gradient
do k=1,nfo
  W(1:npt)=ZERO
  do l=1+Glob_ProcID,nfa,Glob_NumOfProcs
    t=Glob_c(l)
    do m=1,npt
      W(m)=W(m)+t*(Glob_D(m,k,l)-Evalue*Glob_D(m+npt,k,l))
	enddo
  enddo
  t=Glob_c(k+nfru)
  do m=1,npt
    Glob_WkGR((k-1)*npt+m)=2*real(conjg(t)*W(m),dprec)
  enddo
  do m=1+Glob_ProcID,npt,Glob_NumOfProcs
    Glob_WkGR((k-1)*npt+m)=Glob_WkGR((k-1)*npt+m)-t*conjg(t)*(Glob_D(m,k,k+nfru)-Evalue*Glob_D(m+npt,k,k+nfru))
  enddo
enddo

if ((Glob_OverlapPenaltyAllowed).and.(nfa/=1)) then
  call ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty,Glob_WkGR)
  Evalue=Evalue+Glob_TotalOverlapPenalty
endif
call MPI_ALLREDUCE(Glob_WkGR,Gradient,nfo*npt,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

Glob_EnergyGBCounter=Glob_EnergyGBCounter+1

end subroutine EnergyGB



function EnergyIA(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyIA computes the energy of the system under 
!consideration. Routine GSEPIIS from module linalg is used to 
!solve GSEP. The function computes the energy level which is
!the closest one to the value of the global variable 
!Glob_ApproxEnergy. The calculations are done with a basis 
!of Nmax functions. It is assumed that nonlinear parameters
!of the basis functions are stored in global array 
!Glob_NonlinParam and all matrix elements of the
!first Nmin-1 functions have already been calculated and stored in 
!proper global arrays. If no matrix elements have been calculated,
!one should set Nmin=0. In the case when all matrix elements have 
!been computed and only the solution of GSEP is needed, one 
!should set AreMatElemNeeded=.false.
!Upon successfull exit ErrorCode should be equal to 0. In the
!case of nonzero ErrorCode please see the description of GSEPIIS
!as ErrorCode is simply passed from that routine to EnergyIA

real(dprec)    EnergyIA
!Arguments:
integer        Nmin,Nmax
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables:
real(dprec)    Evalue
integer        NumOfIterations


if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
if (Nmax==1) then 
  EnergyIA=real(Glob_H(1,1))+Glob_ApproxEnergy
  NumOfIterations=1
  ErrorCode=0
else
  Glob_c(1:Nmax)=Glob_LastEigvector(1:Nmax)
  call GHEPIIS(Nmin,Nmax,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim,Glob_ApproxEnergy,Glob_c,Glob_WorkForGHEPIIS,Glob_EigvalTol, &
           Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGHEPIIS,-1,NumOfIterations,ErrorCode)
    !GHEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
    !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
  if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol       
  if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol    
  if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_OverlapPenaltyAllowed) then
    EnergyIA=Evalue+Glob_TotalOverlapPenalty   
  else  
    EnergyIA=Evalue
  endif 
endif

Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
Glob_EnergyIACounter=Glob_EnergyIACounter+1

end function EnergyIA



function EnergyIAM(Nmin,Nmax,AreMatElemNeeded,ErrorCode)  
!Function EnergyIAM is pretty much the same function as 
!EnergyIAM. The only exception is that it also generates
!a properly normalized eigenvector, which gives correct 
!linear coefficients (stored in Glob_c).

real(dprec)    EnergyIAM
!Arguments:
integer        Nmin,Nmax
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables:
real(dprec)    Evalue
integer        NumOfIterations


if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
if (Nmax==1) then 
  EnergyIAM=real(Glob_H(1,1))+Glob_ApproxEnergy
  Glob_c(1)=ONE
  NumOfIterations=1
  ErrorCode=0
else
  Glob_c(1:Nmax)=Glob_LastEigvector(1:Nmax)
  call GHEPIIS(Nmin,Nmax,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim,Glob_ApproxEnergy,Glob_c,Glob_WorkForGHEPIIS,Glob_EigvalTol, &
           Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGHEPIIS,0,NumOfIterations,ErrorCode)
    !GHEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
    !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
  if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol       
  if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol   
  if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_OverlapPenaltyAllowed) then
    EnergyIAM=Evalue+Glob_TotalOverlapPenalty   
  else  
    EnergyIAM=Evalue
  endif 
endif

Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
Glob_EnergyIACounter=Glob_EnergyIACounter+1

end function EnergyIAM  



subroutine EnergyIB(Evalue,Gradient,AreMatElemNeeded,ErrorCode)
!subroutine EnergyIB computes the energy value and the
!gradient with respect to the nonlinear parameters of the last 
!Glob_nfo functions of the system under consideration.
!Routine GSEPIIS from module linalg is used to 
!solve GSEP. The subroutine computes the energy level which is
!the closest one to the value of the global variable 
!Glob_ApproxEnergy, as well as the corresponding gradient. 
!It is assumed that nonlinear parameters of the basis functions are
!stored in global array Glob_NonlinParam and all the matrix
!elements of the first Glob_nfru functions have been computed
!and stored in proper global arrays. If no matrix elements have
!been calculated, one should first call ComputeMatElemAndDeriv
!and calculate the matrix elements of the first Glob_nfru
!functions. In case when all matrix elements have been computed
!and it is only necessary to solve GSEP and find gradient, one
!needs to set AreMatElemNeeded=.false.
!Upon successfull exit ErrorCode should be equal to 0. In the
!case of nonzero ErrorCode please see the description of GSEPIIS
!as ErrorCode is simply passed from that routine to EnergyIA.
!The data in Gradient is ordered in the following manner:
!Gradient=(dEdvechLi,dEdvechBi,dEdvechL{i+1},...,dEdvechLj,dEdvechBj)
!where i=Glob_nfru+1, j=Glob_nfa

!Arguments:
real(dprec)    Evalue
real(dprec)    Gradient(Glob_npt*Glob_nfo)
logical        AreMatElemNeeded
integer        ErrorCode

!Local variables:
integer        nfo,nfa,nfru,npt
integer        i,j,k,l,m,nbands,leftover
logical        oddband
integer        NumOfIterations
complex(dprec)    W(Glob_npt_MaxAllowed),t
real(dprec)    pen_coeff

nfo=Glob_nfo
nfa=Glob_nfa
npt=Glob_npt
nfru=Glob_nfru

if (AreMatElemNeeded) call ComputeMatElemAndDeriv(nfru+1,nfa)

if (nfa==1) then
  Evalue=real(Glob_H(1,1))+Glob_ApproxEnergy
  Glob_c(1)=ONE
  NumOfIterations=1
  ErrorCode=0
else
  Glob_c(1:nfa)=Glob_LastEigvector(1:nfa)
  call GHEPIIS(nfru+1,nfa,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim,Glob_ApproxEnergy,Glob_c,Glob_WorkForGHEPIIS,Glob_EigvalTol, &
           Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGHEPIIS,0,NumOfIterations,ErrorCode)
    !GHEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
    !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
  Glob_c(1:nfa)=Glob_LastEigvector(1:nfa)  
  if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol       
  if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol    
  call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

!Computing gradient
do k=1,nfo
  W(1:npt)=ZERO
  do l=1+Glob_ProcID,nfa,Glob_NumOfProcs
    t=Glob_c(l)
    do m=1,npt
      W(m)=W(m)+t*(Glob_D(m,k,l)-Evalue*Glob_D(m+npt,k,l))
	enddo
  enddo
  t=Glob_c(k+nfru)
  do m=1,npt
    Glob_WkGR((k-1)*npt+m)=2*real(conjg(t)*W(m),dprec)
  enddo
  do m=1+Glob_ProcID,npt,Glob_NumOfProcs
    Glob_WkGR((k-1)*npt+m)=Glob_WkGR((k-1)*npt+m)-t*conjg(t)*(Glob_D(m,k,k+nfru)-Evalue*Glob_D(m+npt,k,k+nfru))
  enddo
enddo

if ((Glob_OverlapPenaltyAllowed).and.(nfa/=1)) then
  call ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty,Glob_WkGR)
  Evalue=Evalue+Glob_TotalOverlapPenalty
endif
call MPI_ALLREDUCE(Glob_WkGR,Gradient,nfo*npt,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
Glob_EnergyIBCounter=Glob_EnergyIBCounter+1

end subroutine EnergyIB



subroutine ReadSwapFileAndDistributeData(IsSwapFileOK)
!This subroutine reads data from a swap file (if it is 
!available, and if global variable Glob_UseSwapFile=.true.)
!and sends the read data to all processes. If the result is
!success then the value of a logical variable IsSwapFileOK
!is .true. on exit. It is assumed that the current basis
!size is equal to Glob_CurrBasisSize.

!Arguments:
logical IsSwapFileOK   
!Local variables:
integer i,j,OpenFileErr

IsSwapFileOK=.false.
if (Glob_UseSwapFile) then
  if (Glob_ProcID==0) then
    open(1,file=Glob_SwapFileName,form='binary',status='old',iostat=OpenFileErr)
	if (OpenFileErr==0) then
      read(1) i
	  if (i==Glob_CurrBasisSize) then
        !Reading H and S from a matrix stored in file
		!Remember that the lower part (including the diagonal)
		!contains elements of H, while the upper part contains S
        write(*,'(1x,a33)',advance='no') 'Reading H and S from swap file...'
		do j=1,Glob_CurrBasisSize
          if (OpenFileErr==0) then
             read(1,iostat=OpenFileErr) Glob_H(1:Glob_CurrBasisSize,j)
		  endif
		enddo  
        !reading diagS
        if (OpenFileErr==0) read(1,iostat=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize) 
		if (OpenFileErr==0) then
		  IsSwapFileOK=.true.
	      write(*,*) 'completed'
		  close(1)
        else
          write(*,*) 'failed'
		endif
		!erase information in swap file to free disc space
        open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
		write(1) 'Swap file is empty'
		close(1)
      endif
	endif
  endif
endif 
call MPI_BCAST(IsSwapFileOK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!If swap file is OK then send the data to all processes
if (IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a35)',advance='no') 'Sending H and S to all processes...' 
  call MPI_BCAST(Glob_H,Glob_HSLeadDim*Glob_HSLeadDim,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(Glob_diagS,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  !Remember that the lower part (including the diagonal)
  !contains elements of H, while the upper part contains S 
  
  !Restoring proper storage of data for Glob_GSEPSolutionMethod='G' 
  if (Glob_GSEPSolutionMethod=='G') then
    do i=1,Glob_CurrBasisSize
      do j=1,i-1
        Glob_S(i,j)=Glob_H(j,i)
	  enddo
	  Glob_diagH(i)=real(Glob_H(i,i),dprec)
    enddo
  endif
  !Restoring proper storage of data for Glob_GSEPSolutionMethod='I'
  if (Glob_GSEPSolutionMethod=='I') then
    do i=1,Glob_CurrBasisSize
      do j=1,i-1
        Glob_S(i,j)=Glob_H(j,i)
        Glob_S(j,i)=conjg(Glob_H(j,i))
	  enddo
	  Glob_S(i,i)=ONE
    enddo
    do i=1,Glob_CurrBasisSize
      Glob_H(i,i)=Glob_H(i,i)-Glob_ApproxEnergy
      do j=i+1,Glob_CurrBasisSize
        Glob_H(j,i)=Glob_H(j,i)-Glob_ApproxEnergy*Glob_S(j,i)
      enddo
    enddo
  endif    
    
  if (Glob_ProcID==0) write(*,*) 'completed' 
else
  if (Glob_ProcID==0) then 
    write(*,*) 'Matrices H and S were not read from swap file'
    write(*,*) 'All H and S matrix elements need be (re)computed'
  endif
endif

end subroutine ReadSwapFileAndDistributeData


subroutine StoreMatricesInSwapFile()
!This subroutine stores H and S matrices, including
!their diagonals, in a swap file. This is done only
!if the value of global variable Glob_UseSwapFile is .true.

!Local variables:
integer i,j,OpenFileErr

if (Glob_UseSwapFile) then
!If swap file is allowed to use then write H and S matrix
!elements into it
  if (Glob_CurrBBOPStep/=Glob_NumOfBBOPSteps) then
    if (Glob_ProcID==0) then
      open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
	  if (OpenFileErr==0) then
        write(1) Glob_CurrBasisSize
        !We store a matrix whose lower part (including the diagonal)
	    !contains elements of H, while the upper part contains S.
	    !Also, we store the diagonal of S.

        if (Glob_GSEPSolutionMethod=='G') then
          do i=1,Glob_CurrBasisSize
            do j=1,i-1
              Glob_H(j,i)=Glob_S(i,j)
	        enddo
	        Glob_H(i,i)=Glob_diagH(i)
          enddo
        endif
        if (Glob_GSEPSolutionMethod=='I') then
          do i=1,Glob_CurrBasisSize
            do j=1,i-1
              Glob_H(j,i)=Glob_S(i,j)
	        enddo
	        Glob_H(i,i)=Glob_H(i,i)+Glob_ApproxEnergy
	        do j=i+1,Glob_CurrBasisSize
	          Glob_H(j,i)=Glob_H(j,i)+Glob_ApproxEnergy*Glob_S(j,i)
	        enddo
          enddo       
        endif  
                
        write(*,'(1x,a33)',advance='no') 'Writing H and S into swap file... '
		do i=1,Glob_CurrBasisSize
          if (OpenFileErr==0) then
            write(1,iostat=OpenFileErr) Glob_H(1:Glob_CurrBasisSize,i)
		  endif
		enddo         
        if (OpenFileErr==0) then
	      write(1,iostat=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize) 
        endif
        if (OpenFileErr==0) then
	      write(*,*) 'completed'
	    else
          write(*,*) 'failed to complete'
	    endif
      endif
	  close(1)
    endif
  endif
endif  

end subroutine StoreMatricesInSwapFile



subroutine ReadHessianFile(V,IVLMAT,D,nvar,FileName,IsHessFileOK)
!This subroutine reads Hessian from file FileName (if it exists
!and contains appropriate amount of data). nvar stands for the 
!total number of variables in the problem. The Hessian is 
!placed in array V, starting from element V(IVLMAT).
!In case when the global parameter Glob_FullOptSaveD=.true. 
!the subroutine also reads scaling vector D from the same 
!file. It assumed that this scaling vector is placed after 
!Hessian data in the file.
!If the result is a success then the value of a logical variable 
!IsHessFileOK is .true. on exit.

!Arguments:
real(dprec) V(*),D(*)
integer     IVLMAT,nvar
character(Glob_FileNameLength) FileName
logical     IsHessFileOK
!Local variables:
integer i,j,OpenFileErr

IsHessFileOK=.false.
if (Glob_ProcID==0) then
  open(1,file=FileName,form='binary',status='old',iostat=OpenFileErr)
  if (OpenFileErr==0) then
    read(1,iostat=OpenFileErr) j
    if ((j==nvar).and.(OpenFileErr==0)) then
      write(*,'(1x,a28)',advance='no') 'Reading Hessian from file...'
      read(1,iostat=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)
	  if ((OpenFileErr==0).and.(Glob_FullOptSaveD)) read(1,iostat=OpenFileErr) D(1:nvar)
	  if (OpenFileErr==0) then
		IsHessFileOK=.true.
	    write(*,*) 'done'
		close(1)
      else
        write(*,*) 'failed'
        write(*,*) 'Warning: Default Hessian initialization must be used'
	  endif
	else
      write(*,*) 'Warning: Hessian file ',FileName
	  write(*,*) 'is inconsistent with the dimension of the current optimization problem' 
	  write(*,*) 'and will not be read'
	  write(*,*) 'Default Hessian initialization must be used'      
	endif
  else
    write(*,*) 'Warning: Cannot open Hessian file ',FileName
	write(*,*) 'Default Hessian initialization must be used'
  endif
endif

end subroutine ReadHessianFile



subroutine SaveHessianFile(V,IVLMAT,D,nvar,FileName,IsSuccess)
!This subroutine saves Hessian to file FileName. nvar stands for the 
!total number of variables in the problem. The Hessian is taken from
!V, starting from element V(IVLMAT).
!In case when the golbal parameter Glob_FullOptSaveD=.true. the 
!subroutine also saves the scaling vector D (used by the optimization 
!routine). Vector D is placed after V, i.e. in the end of the file, so
!that iven if user happens to continue the same optimization problem with 
!Glob_FullOptSaveD changed from .true. to .false., the routine
!will run without failure.
!If the file is successfully saved then the value of a logical variable 
!IsSuccess is .true. on exit.

!Arguments:
real(dprec) V(*),D(*)
integer     IVLMAT,nvar
character(Glob_FileNameLength) FileName
logical     IsSuccess
!Local variables:
integer i,OpenFileErr

IsSuccess=.false.
if (Glob_ProcID==0) then
  open(1,file=FileName,form='binary',status='replace',iostat=OpenFileErr)
  if (OpenFileErr==0) then
    write(1,iostat=OpenFileErr) nvar
    if (OpenFileErr==0) then
      write(*,'(1x,a17)',advance='no') 'Saving Hessian...'
      write(1,iostat=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)
	  if (Glob_FullOptSaveD) write(1,iostat=OpenFileErr) D(1:nvar)
	  close(1)
	  if (OpenFileErr==0) then
		IsSuccess=.true.
	    write(*,*) 'done'
      else
        write(*,*) 'failed'
	  endif     
	endif
  endif
  if (.not.IsSuccess) write(*,*) 'Warning: Hessian was not saved'
endif

end subroutine SaveHessianFile



subroutine PermuteFunctions(fb,fe,FuncNumTemp,NonlinParamTemp)
!Subroutine PermuteFunctions permutes basis functions in such 
!a way that functions fb through fe go to the very end of the
!basis (which makes their optimization more time efficient). 
!The subroutine permutes the nonlinear parameters, Glob_NonlinParam, 
!and function numbers, Glob_FuncNum.
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp and NonlinParamTemp. The length of FuncNumTemp
!should be at least fe-fb+1, while the lenght of NonlinParamTemp must be
!at least (fe-fb+1)*Glob_npt

!Arguments:
integer      fb,fe
integer      FuncNumTemp(*)
real(dprec)      NonlinParamTemp(Glob_npt,*)
!Local variables:
integer  i,fbn,nfco

nfco=fe-fb+1
fbn=fb+Glob_CurrBasisSize-fe

NonlinParamTemp(1:Glob_npt,1:nfco)=Glob_NonlinParam(1:Glob_npt,fb:fe)
do i=fb,fbn-1
  Glob_NonlinParam(1:Glob_npt,i)=Glob_NonlinParam(1:Glob_npt,nfco+i)
enddo
Glob_NonlinParam(1:Glob_npt,fbn:Glob_CurrBasisSize)=NonlinParamTemp(1:Glob_npt,1:nfco)

FuncNumTemp(1:nfco)=Glob_FuncNum(fb:fe)
do i=fb,fbn-1
  Glob_FuncNum(i)=Glob_FuncNum(nfco+i)
enddo
Glob_FuncNum(fbn:Glob_CurrBasisSize)=FuncNumTemp(1:nfco)

end subroutine PermuteFunctions



subroutine PermuteFunctions2(fb1,fe1,fe2,FuncNumTemp,NonlinParamTemp)
!Subroutine PermuteFunctions2 permutes basis functions in such 
!a way that the set of functions fb1 through fe1 exchanges its positions
!with the set of function fb2 through fe2, where fb2=fe1+1. 
!The subroutine permutes the nonlinear parameters, Glob_NonlinParam, 
!and function numbers, Glob_FuncNum.
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp and NonlinParamTemp. The length of FuncNumTemp
!should be at least fe1-fb1+1, while the lenght of NonlinParamTemp must be
!at least (fe1-fb1+1)*Glob_npt

!Arguments:
integer      fb1,fe1,fe2
integer      FuncNumTemp(*)
real(dprec)      NonlinParamTemp(Glob_npt,*)
!Local variables:
integer  i,fbn,k

k=fe1-fb1+1
fbn=fb1+fe2-fe1

NonlinParamTemp(1:Glob_npt,1:k)=Glob_NonlinParam(1:Glob_npt,fb1:fe1)
do i=fb1,fbn-1
  Glob_NonlinParam(1:Glob_npt,i)=Glob_NonlinParam(1:Glob_npt,k+i)
enddo
Glob_NonlinParam(1:Glob_npt,fbn:fe2)=NonlinParamTemp(1:Glob_npt,1:k)
  
FuncNumTemp(1:k)=Glob_FuncNum(fb1:fe1)
do i=fb1,fbn-1
  Glob_FuncNum(i)=Glob_FuncNum(k+i)
enddo
Glob_FuncNum(fbn:fe2)=FuncNumTemp(1:k)

end subroutine PermuteFunctions2



subroutine PermuteMatrixElements(fb,fe,TempR,TempC)
!This subroutine permutes matrix elements of H and S in such a 
!way that functions fb through fe go to the very end of the basis 
!(which makes their optimization more time efficient). The subroutine 
!permutes elements of arrays Glob_H, Glob_S, Glob_diagH, Glob_diagS
!This subroutine requires workspace, which must be provided by
!arrays TempR and TempC. The length of these araays should be at 
!least fe-fb+1. 
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper 
!triangles (excluding the diagonal) are used as workspace.   

!Arguments:
integer       fb,fe
real(dprec)       TempR(*)
complex(dprec)    TempC(*)
!Local variables:
integer  i,j,fbn,nfco,fep,k,q,p,r,s

nfco=fe-fb+1
fbn=fb+Glob_CurrBasisSize-fe

!First we do the diagonals of H and S
if (Glob_GSEPSolutionMethod=='G') then
  TempR(1:nfco)=Glob_diagH(fb:fe)
  do i=fb,fbn-1
    Glob_diagH(i)=Glob_diagH(nfco+i)
  enddo
  Glob_diagH(fbn:Glob_CurrBasisSize)=TempR(1:nfco)
else
  do i=fb,fe
    TempC(i-fb+1)=Glob_H(i,i)
  enddo
  do i=fb,fbn-1
    Glob_H(i,i)=Glob_H(nfco+i,nfco+i)
  enddo
  do i=fbn,Glob_CurrBasisSize
    Glob_H(i,i)=TempC(i-fbn+1)
  enddo
endif 
    
TempR(1:nfco)=Glob_diagS(fb:fe)
do i=fb,fbn-1
  Glob_diagS(i)=Glob_diagS(nfco+i)
enddo
Glob_diagS(fbn:Glob_CurrBasisSize)=TempR(1:nfco)

!Now we permute off-diagonal elements
fep=fe+1
k=Glob_CurrBasisSize-fe
q=Glob_CurrBasisSize+fep
p=fb+Glob_CurrBasisSize-fe-1
s=Glob_CurrBasisSize+fb
r=s-1

do i=1,fb-1
  TempC(1:nfco)=Glob_H(fb:fe,i)
  do j=fb,fbn-1
    Glob_H(j,i)=Glob_H(j+nfco,i)
  enddo
  Glob_H(fbn:Glob_CurrBasisSize,i)=TempC(1:nfco)
enddo
do i=fep,Glob_CurrBasisSize-1
  Glob_H(fep:q-i-1,q-i)=Glob_H(i+1:Glob_CurrBasisSize,i)
enddo
do i=fe+1,Glob_CurrBasisSize
  do j=fb,fe
    Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
enddo
do i=fb,fe-1
  Glob_H(i+k+1:Glob_CurrBasisSize,i+k)=Glob_H(i+1:fe,i) 
enddo
Glob_H(fbn:Glob_CurrBasisSize,fb:fb+k-1)=Glob_H(fb:fe,fe+1:Glob_CurrBasisSize)
do i=fb,p-1
  Glob_H(i+1:p,i)=Glob_H(fep:r-i,s-i)
enddo

do i=1,fb-1
  TempC(1:nfco)=Glob_S(fb:fe,i)
  do j=fb,fbn-1
    Glob_S(j,i)=Glob_S(j+nfco,i)
  enddo
  Glob_S(fbn:Glob_CurrBasisSize,i)=TempC(1:nfco)
enddo
do i=fep,Glob_CurrBasisSize-1
  Glob_S(fep:q-i-1,q-i)=Glob_S(i+1:Glob_CurrBasisSize,i)
enddo
do i=fep,Glob_CurrBasisSize
  do j=fb,fe
    Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
enddo
do i=fb,fe-1
  Glob_S(i+k+1:Glob_CurrBasisSize,i+k)=Glob_S(i+1:fe,i) 
enddo
Glob_S(fbn:Glob_CurrBasisSize,fb:fb+k-1)=Glob_S(fb:fe,fe+1:Glob_CurrBasisSize)
do i=fb,p-1
  Glob_S(i+1:p,i)=Glob_S(fep:r-i,s-i)
enddo
if (Glob_GSEPSolutionMethod=='I') then
  !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
  !upper triangle of Glob_S
  do i=fb,Glob_CurrBasisSize
    do j=1,i-1
      Glob_S(j,i)=conjg(Glob_S(i,j))
    enddo
  enddo
endif

end subroutine PermuteMatrixElements



subroutine PermuteMatrixElements2(fb1,fe1,fe2,TempR,TempC)
!This subroutine permutes matrix elements of H and S in such a 
!way that the set of functions fb1 through fe1 exchanges its 
!positions with the set of function fb2 through fe2, where 
!fb2=fe1+1.  
!The subroutine permutes elements of arrays Glob_H, Glob_S, 
!Glob_diagH, Glob_diagS. This subroutine requires workspace, 
!which must be provided by arrays TempR and TempC. The length of 
!these araays should be at least fe1-fb1+1. 
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper 
!triangles (excluding the diagonal) are used as workspace. 
!Arguments:
integer       fb1,fe1,fe2
real(dprec)       TempR(*)
complex(dprec)    TempC(*)
!Local variables:
integer  i,j,fn,t,fe1p,k,q,p,r,s,fb2

t=fe1-fb1+1
fn=fb1+fe2-fe1
fb2=fe1+1

!First we do the diagonals of H and S
if (Glob_GSEPSolutionMethod=='G') then
  TempR(1:t)=Glob_diagH(fb1:fe1)
  do i=fb1,fn-1
    Glob_diagH(i)=Glob_diagH(t+i)
  enddo
  Glob_diagH(fn:fe2)=TempR(1:t)
else
  do i=fb1,fe1
    TempC(i-fb1+1)=Glob_H(i,i)
  enddo
  do i=fb1,fn-1
    Glob_H(i,i)=Glob_H(t+i,t+i)
  enddo
  do i=fn,fe2
    Glob_H(i,i)=TempC(i-fn+1)
  enddo
endif
    
TempR(1:t)=Glob_diagS(fb1:fe1)
do i=fb1,fn-1
  Glob_diagS(i)=Glob_diagS(t+i)
enddo
Glob_diagS(fn:fe2)=TempR(1:t)

!Now we permute off-diagonal elements
fe1p=fe1+1
k=fe2-fe1
q=fe2+fe1p
p=fb1+fe2-fe1-1
s=fe2+fb1
r=s-1

do i=1,fb1-1
  TempC(1:t)=Glob_H(fb1:fe1,i)
  do j=fb1,fn-1
    Glob_H(j,i)=Glob_H(t+j,i)
  enddo
  Glob_H(fn:fe2,i)=TempC(1:t)
enddo
do i=fe1p,fe2-1
  Glob_H(fe1p:q-i-1,q-i)=Glob_H(i+1:fe2,i)
enddo
do i=fe1+1,fe2
  do j=fb1,fe1
    Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
enddo
do i=fb1,fe1-1
  Glob_H(i+k+1:fe2,i+k)=Glob_H(i+1:fe1,i) 
enddo
Glob_H(fn:fe2,fb1:fb1+k-1)=Glob_H(fb1:fe1,fe1+1:fe2)
do i=fb1,p-1
  Glob_H(i+1:p,i)=Glob_H(fe1p:r-i,s-i)
enddo
do i=fe2+1,Glob_CurrBasisSize
  TempC(1:t)=Glob_H(i,fb1:fe1)
  do j=fb1,fn-1
    Glob_H(i,j)=Glob_H(i,t+j)
  enddo
  Glob_H(i,fn:fe2)=TempC(1:t)
enddo

do i=1,fb1-1
  TempC(1:t)=Glob_S(fb1:fe1,i)
  do j=fb1,fn-1
    Glob_S(j,i)=Glob_S(t+j,i)
  enddo
  Glob_S(fn:fe2,i)=TempC(1:t)
enddo
do i=fe1p,fe2-1
  Glob_S(fe1p:q-i-1,q-i)=Glob_S(i+1:fe2,i)
enddo
do i=fe1+1,fe2
  do j=fb1,fe1
    Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
enddo
do i=fb1,fe1-1
  Glob_S(i+k+1:fe2,i+k)=Glob_S(i+1:fe1,i) 
enddo
Glob_S(fn:fe2,fb1:fb1+k-1)=Glob_S(fb1:fe1,fe1+1:fe2)
do i=fb1,p-1
  Glob_S(i+1:p,i)=Glob_S(fe1p:r-i,s-i)
enddo
do i=fe2+1,Glob_CurrBasisSize
  TempC(1:t)=Glob_S(i,fb1:fe1)
  do j=fb1,fn-1
    Glob_S(i,j)=Glob_S(i,t+j)
  enddo
  Glob_S(i,fn:fe2)=TempC(1:t)
enddo
if (Glob_GSEPSolutionMethod=='I') then
  !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
  !upper triangle of Glob_S
  do i=fb1,Glob_CurrBasisSize
    do j=1,i-1
      Glob_S(j,i)=conjg(Glob_S(i,j))
    enddo
  enddo
endif

end subroutine PermuteMatrixElements2



subroutine ReverseFuncOrder(fb,fe)
!Subroutine ReverseFuncOrder changes the order of basis functions 
!fb through fe to reverse. 
!The subroutine permutes the nonlinear parameters, Glob_NonlinParam, 
!and function numbers, Glob_FuncNum.

!Arguments:
integer  fb,fe
!Local variables:
real(dprec) temp(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
integer i,j,f,t,fbm,fep

f=(fe-fb+1)/2 !integer division!
fbm=fb-1
fep=fe+1
do i=1,f
  temp(1:Glob_npt)=Glob_NonlinParam(1:Glob_npt,fbm+i)
  Glob_NonlinParam(1:Glob_npt,fbm+i)=Glob_NonlinParam(1:Glob_npt,fep-i)
  Glob_NonlinParam(1:Glob_npt,fep-i)=temp(1:Glob_npt)
enddo
do i=1,f
  t=Glob_FuncNum(fbm+i)
  Glob_FuncNum(fbm+i)=Glob_FuncNum(fep-i)
  Glob_FuncNum(fep-i)=t
enddo

end subroutine ReverseFuncOrder



subroutine ReverseMatElemOrder(fb,fe)
!Subroutine ReverseMatElemOrder changes the order of matrix elements
!corresponding to basis functions fb through fe to reverse. 
!The subroutine permutes elements of arrays Glob_H, Glob_S, Glob_diagH, 
!Glob_diagS. Note that it uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonal). Upper triangles are not 
!referenced.

!Arguments:
integer  fb,fe
!Local variables:
integer i,j,f,fbm,fep,ff
real(dprec)     r
complex(dprec)  c

f=(fe-fb+1)/2 !integer division!
fbm=fb-1
fep=fe+1
ff=fe+fb
!Diagonal elements
if (Glob_GSEPSolutionMethod=='G') then
  do i=1,f
    r=Glob_diagH(fbm+i)
    Glob_diagH(fbm+i)=Glob_diagH(fep-i)
    Glob_diagH(fep-i)=r
  enddo
else
  do i=1,f
    c=Glob_H(fbm+i,fbm+i)
    Glob_H(fbm+i,fbm+i)=Glob_H(fep-i,fep-i)
    Glob_H(fep-i,fep-i)=c
  enddo
endif
do i=1,f
  r=Glob_diagS(fbm+i)
  Glob_diagS(fbm+i)=Glob_diagS(fep-i)
  Glob_diagS(fep-i)=r
enddo
!Off-diagonal elements
do i=1,fbm
  do j=1,f
    c=Glob_H(fbm+j,i)
    Glob_H(fbm+j,i)=Glob_H(fep-j,i)
    Glob_H(fep-j,i)=c
  enddo
enddo
do i=1,f
  do j=fb+i,fe-i
    c=Glob_H(j,fbm+i)
    Glob_H(j,fbm+i)=conjg(Glob_H(fep-i,ff-j))
	Glob_H(fep-i,ff-j)=conjg(c)
  enddo
  Glob_H(j,fbm+i)=conjg(Glob_H(j,fbm+i))
enddo
do i=fb,fb+f-1
  do j=fep,Glob_CurrBasisSize
    c=Glob_H(j,i)
    Glob_H(j,i)=Glob_H(j,ff-i)
    Glob_H(j,ff-i)=c
  enddo
enddo
do i=1,fbm
  do j=1,f
    c=Glob_S(fbm+j,i)
    Glob_S(fbm+j,i)=Glob_S(fep-j,i)
    Glob_S(fep-j,i)=c
  enddo
enddo
do i=1,f
  do j=fb+i,fe-i
    c=Glob_S(j,fbm+i)
    Glob_S(j,fbm+i)=conjg(Glob_S(fep-i,ff-j))
	Glob_S(fep-i,ff-j)=conjg(c)
  enddo
  Glob_S(j,fbm+i)=conjg(Glob_S(j,fbm+i))
enddo
do i=fb,fb+f-1
  do j=fep,Glob_CurrBasisSize
    c=Glob_S(j,i)
    Glob_S(j,i)=Glob_S(j,ff-i)
    Glob_S(j,ff-i)=c
  enddo
enddo
if (Glob_GSEPSolutionMethod=='I') then
  !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
  !upper triangle of Glob_S
  do i=fb,fe
    do j=1,i-1
      Glob_S(j,i)=conjg(Glob_S(i,j))
    enddo
  enddo
endif

end subroutine ReverseMatElemOrder


subroutine SortBasisFuncAndMatElem(fb,fe,FuncNumTemp,NonlinParamTemp,TempR,TempC)
!Subroutine SortBasisFuncAndMatElem permutes basis functions
!fb through fe so that they are sorted in decreasing order
!(that is the values of Glob_FuncNum(i) decrease). It also permutes
!the corresponding matrix elements.
!The subroutine permutes the nonlinear parameters, Glob_NonlinParam, 
!function numbers, Glob_FuncNum, the elements of arrays Glob_H, Glob_S,
!Glob_diagH, glob_diagS. It is assumed that the set of
!values Glob_FuncNum(fb:fe) ranges from some minimal value
!to the maximal one with no gaps (for example:  5,8,6,7,9). If
!this condition is not satisfied the subroutine will fail without
!warning.
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp, NonlinParamTemp, TempR, TempC. The length of 
!FuncNumTemp, TempR, TempC should be at least fe-fb+1, while the 
!lenght of NonlinParamTemp must be at least (fe-fb+1)*Glob_npt.
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper 
!triangles (excluding the diagonal) are used as workspace.  

!Arguments:
integer       fb,fe
integer       FuncNumTemp(*)
real(dprec)       NonlinParamTemp(Glob_npt,*)
real(dprec)       TempR(*)
complex(dprec)    TempC(*)
!Local variables:
integer  i,j,k,mf,fep,nfco,fbm,cbs,fi,fj

cbs=Glob_CurrBasisSize
nfco=fe-fb+1
fep=fe+1
fbm=fb-1
mf=minval(Glob_FuncNum(fb:fe))
k=fbm+nfco+mf
!First we sort out nonlinear parameters
do i=fb,fe
  NonlinParamTemp(1:Glob_npt,nfco+mf-Glob_FuncNum(i))=Glob_NonlinParam(1:Glob_npt,i)
enddo
Glob_NonlinParam(1:Glob_npt,fb:fe)=NonlinParamTemp(1:Glob_npt,1:nfco)
!Then we sort out the diagonal matrix elements
if (Glob_GSEPSolutionMethod=='G') then
  do i=fb,fe
    TempR(nfco+mf-Glob_FuncNum(i))=Glob_diagH(i)
  enddo
  Glob_diagH(fb:fe)=TempR(1:nfco)
else
  do i=fb,fe
    TempC(nfco+mf-Glob_FuncNum(i))=Glob_H(i,i)
  enddo
  do i=fb,fe
    Glob_H(i,i)=TempC(i-fb+1) 
  enddo   
endif
do i=fb,fe
  TempR(nfco+mf-Glob_FuncNum(i))=Glob_diagS(i)
enddo
Glob_diagS(fb:fe)=TempR(1:nfco)
!Sorting out off-diagonal matrix elements
do i=fb,fe
  Glob_H(1:fbm,k-Glob_FuncNum(i))=Glob_H(i,1:fbm)
enddo
do i=fb,fe
  Glob_H(i,1:fbm)=Glob_H(1:fbm,i)
enddo
do i=fb,fe
  Glob_H(k-Glob_FuncNum(i),fep:cbs)=Glob_H(fep:cbs,i)
enddo
do i=fb,fe
  Glob_H(fep:cbs,i)=Glob_H(i,fep:cbs)
enddo
do i=fb,fe-1
  do j=i+1,fe
    fi=Glob_FuncNum(i)
	fj=Glob_FuncNum(j)
    if (fi<fj) then
       Glob_H(k-fj,k-fi)=conjg(Glob_H(j,i))
	else
       Glob_H(k-fi,k-fj)=Glob_H(j,i)
	endif
  enddo
enddo
do i=fb,fe-1
    Glob_H(i+1:fe,i)=Glob_H(i,i+1:fe)
enddo
do i=fb,fe
  Glob_S(1:fbm,k-Glob_FuncNum(i))=Glob_S(i,1:fbm)
enddo
do i=fb,fe
  Glob_S(i,1:fbm)=Glob_S(1:fbm,i)
enddo
do i=fb,fe
  Glob_S(k-Glob_FuncNum(i),fep:cbs)=Glob_S(fep:cbs,i)
enddo
do i=fb,fe
  Glob_S(fep:cbs,i)=Glob_S(i,fep:cbs)
enddo
do i=fb,fe-1
  do j=i+1,fe
    fi=Glob_FuncNum(i)
	fj=Glob_FuncNum(j)
    if (fi<fj) then
       Glob_S(k-fj,k-fi)=conjg(Glob_S(j,i))
	else
       Glob_S(k-fi,k-fj)=Glob_S(j,i)
	endif
  enddo
enddo
do i=fb,fe-1
    Glob_S(i+1:fe,i)=Glob_S(i,i+1:fe)
enddo
if (Glob_GSEPSolutionMethod=='I') then
  !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
  !upper triangle of Glob_S
  do i=fb,fe
    do j=1,i-1
      Glob_S(j,i)=conjg(Glob_S(i,j))
    enddo
  enddo
endif
!At last we change the order of basis functions
do i=fb,fe
   Glob_FuncNum(i)=mf+fe-i
enddo

end subroutine SortBasisFuncAndMatElem



subroutine GetOverlapStatistics(Nmin,Nmax,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
!Subroutine GetOverlapStatistics determines the largest by magnitude overlap, the smallest
!by magnitude overalap, and the average overlap magnitude for basis functions 
!ranging from Nmin to Nmax. 
!Arguments:
integer         Nmin,Nmax
complex(dprec)  MaxAbsOverlap,MinAbsOverlap
real(dprec)     AverageAbsOverlap
!Local variables
integer      i,j,k,nbands,leftover
logical      oddband
real(dprec)  absMaxAbsOverlap,absMinAbsOverlap,absSji,t

MaxAbsOverlap=ZERO
absMaxAbsOverlap=ZERO
MinAbsOverlap=huge(absMinAbsOverlap)
absMinAbsOverlap=MinAbsOverlap
t=ZERO

do i=1,Nmin-1
  do j=Nmin,Nmax
    absSji=abs(Glob_S(j,i))
    if (absSji>absMaxAbsOverlap) then
       absMaxAbsOverlap=absSji
       MaxAbsOverlap=Glob_S(j,i)
    endif   
    if (absSji<absMinAbsOverlap) then
       absMinAbsOverlap=absSji
       MinAbsOverlap=Glob_S(j,i)
    endif        
    t=t+absSji
  enddo
enddo  
do i=Nmin,Nmax
  do j=i+1,Nmax
    absSji=abs(Glob_S(j,i))
    if (absSji>absMaxAbsOverlap) then
       absMaxAbsOverlap=absSji
       MaxAbsOverlap=Glob_S(j,i)
    endif   
    if (absSji<absMinAbsOverlap) then
       absMinAbsOverlap=absSji
       MinAbsOverlap=Glob_S(j,i)
    endif        
    t=t+absSji
  enddo
enddo    
AverageAbsOverlap=2*t/(Nmax*(Nmax-1)-(Nmin-1)*(Nmin-2)) 
 

end subroutine GetOverlapStatistics



function NumOfRowsToPermForUnitShift(j)
!Function NumOfRowsToPerm is used in cyclic optimization 
!procedure to find the number of function that must be 
!permuted. What it does, it returns the following values:
!
!      j  NumOfRowsToPermForUnitShift(j)
!
!      1           1
!      2           2
!      3           1
!      4           4
!      5           1
!      6           2
!      7           1
!      8           8
!      9           1
!     10           2
!     11           1
!     12           4
!     13           1
!     14           2
!     15           1
!     16          16
!     17           1    
!      .           .   
!      .           . 
!      .           .       .

integer NumOfRowsToPermForUnitShift,j,k
k=1
do while (mod(j,k)==0)
  k=k*2
enddo  
NumOfRowsToPermForUnitShift=k/2

end function NumOfRowsToPermForUnitShift



subroutine BasisEnlG(Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval, &
OverlapThreshold,LinCoeffThreshold) 
!Subroutine BasisEnlG enlarges the basis set from initial
!Kstart-1 functions to Kstop functions by means of trials of 
!randomly selected candidates and then optimizing them. At each 
!step it generates a set of Kstep functions for NTrials times
!based on the existing distribution of nonlinear papameters. 
!Only the set that lowers the energy the most is left. Then this
!set is optimized according to the value of OptimizationType.
!MaxEnergyEval is the limit of the energy evaluations for this
!optimization process. 
!Each newly accepted basis function is checked for pair linear
!dependency with other functions in the basis. If the absolute
!value of the overlap is greater than OverlapThreshold then
!such function is rejected and random trial and optimization 
!take place again until either a good function/functions are 
!generated or a certain limit of such repetitions is reached 
!(the limit is given by the value of global constant 
!Glob_BadOverlapOrLinCoeffLim). If the value of OverlapThreshold 
!passed is negative or zero then no such check is performed. 
!A Similar check is performed for all linear parameters. If
!adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions 
!are rejected and the procedure is repeted until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avaid this check it is enough to set LinCoeffThreshold to a value 
!equal or smaller than zero.
!It is assumed that the basis of Kstart-1 functions has already
!been selected and the corresponding nonlinear parameters are
!stored in proper global array. Subroutines EnergyGA, EnergyGB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)     :: Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval
real(dprec),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
integer      i,j,K,AttemptToGetGoodFunc,ii
integer      np,npt,nfo,nfa,nfru,nfrup1,nvmax,nv
integer      OpenFileErr,ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
logical      IsSwapFileOK,IsEnergyImproved,ExitNeeded
logical      IsOverlapBad,IsAnyLinCoeffBad
integer      wbfu_t,wmu_t,wbfu,wmu,rgm1_counter,rgm2_counter,BlockSizeForZHEGVX
real(dprec)  ms1,ms2
real(dprec)  Evalue,E_init,E_best
real(dprec)  t
real(dprec),allocatable,dimension(:,:)   :: ParSet,ParSetBest
real(dprec),allocatable,dimension(:)     :: x,x_best,grad 
!==================================================== 
!These variables are used when a finite difference gradient is computed               
!real(dprec),allocatable,dimension(:)     ::    fx,fgrad
!real(dprec)                                    deltax,Evalue1
!====================================================
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V,V_init
integer,parameter    :: LIV=60
integer                 IV(LIV), IV_init(LIV)
integer                 LV
integer                 ALG
!Allocatable work space
real(dprec),allocatable,dimension(:)            :: WorkBuffReal
integer,allocatable,dimension(:)                :: WorkBuffInt
type(Glob_HistoryStep),allocatable,dimension(:) :: TempHistory
real(dprec),allocatable,dimension(:,:)          :: TempParam
integer,allocatable,dimension(:)                :: TempFunc

if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine BasisEnlG has started'
  write(*,'(1x,a7,i5)') 'Kstart=',Kstart
  write(*,'(1x,a7,i5)') 'Kstop= ',Kstop
  write(*,'(1x,a7,i5)') 'Kstep= ',Kstep
  write(*,'(1x,a17,i4)') 'OptimizationType=',OptimizationType
  write(*,'(1x,a14,i7)') 'MaxEnergyEval=',MaxEnergyEval
endif

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_nfa=Kstart+Kstep
Glob_nfo=Kstep
Glob_HSLeadDim=Kstop
Glob_HSBuffLen=Kstop*Kstep
np=Glob_np
npt=Glob_npt
nfo=Glob_nfo
nfa=Glob_nfa
nvmax=Kstep*Glob_npt

!Reallocate arrays that contain the information
!about basis functions and optimization process.

if (Glob_UseReallocFile) then
  !Temporarily store the information in a file
  if (Glob_ProcID==0) then
    if (Glob_CurrBasisSize>0) then
      write(*,'(1x,a47)',advance='no') 'Reallocating some arrays using external file...'
      open(1,file=Glob_ReallocFileName,form='binary',status='replace')
      write(1) Glob_History(1:Glob_CurrBasisSize)
      write(1) Glob_FuncNum(1:Glob_CurrBasisSize)
      write(1) Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
      close(1)
	endif
  endif
  deallocate(Glob_NonlinParam)
  deallocate(Glob_FuncNum)
  deallocate(Glob_History)
  allocate(Glob_History(Kstop))
  allocate(Glob_FuncNum(Kstop))
  allocate(Glob_NonlinParam(npt,Kstop))
  if (Glob_ProcID==0) then
    if (Glob_CurrBasisSize>0) then
      open(1,file=Glob_ReallocFileName,form='binary',status='old',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        read(1) Glob_History(1:Glob_CurrBasisSize)
        read(1) Glob_FuncNum(1:Glob_CurrBasisSize)
        read(1,iostat=OpenFileErr) Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
      endif
      close(1) 
      call MPI_BCAST(OpenFileErr,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  endif
  if ((OpenFileErr/=0).and.(Glob_CurrBasisSize>0)) then
    if (Glob_ProcID==0) then
	  write(*,*)
	  write(*,*) 'Error in BasisEnlG: cannot read data from ',Glob_ReallocFileName
    endif
    stop
  endif
  if (Glob_ProcID==0) then
    open(1,file=Glob_ReallocFileName,form='binary',status='replace',iostat=OpenFileErr)
    write(1)  'This temporary file is empty'
    close(1)
  endif
  allocate(WorkBuffReal(Glob_CurrBasisSize))
  allocate(WorkBuffInt(Glob_CurrBasisSize))
  do i=1,Glob_CurrBasisSize
    WorkBuffReal(i)=Glob_History(i)%Energy
  enddo
  call MPI_BCAST(WorkBuffReal,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%Energy=WorkBuffReal(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%CyclesDone
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%CyclesDone=WorkBuffInt(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%InitFuncAtLastStep
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%InitFuncAtLastStep=WorkBuffInt(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%NumOfEnergyEvalDuringFullOpt
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%NumOfEnergyEvalDuringFullOpt=WorkBuffInt(i)
  enddo
  deallocate(WorkBuffReal)
  deallocate(WorkBuffInt)
  call MPI_BCAST(Glob_NonlinParam,npt*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_ProcID==0) then
    write(*,*) 'completed'
  endif
else ! if (Glob_UseReallocFile)
  allocate(TempHistory(Glob_CurrBasisSize))
  allocate(TempFunc(Glob_CurrBasisSize))
  allocate(TempParam(npt,Glob_CurrBasisSize))
  TempHistory(1:Glob_CurrBasisSize)=Glob_History(1:Glob_CurrBasisSize)
  TempFunc(1:Glob_CurrBasisSize)=Glob_FuncNum(1:Glob_CurrBasisSize)
  TempParam(1:npt,1:Glob_CurrBasisSize)=Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
  deallocate(Glob_NonlinParam)
  deallocate(Glob_FuncNum)
  deallocate(Glob_History)
  allocate(Glob_History(Kstop))
  allocate(Glob_FuncNum(Kstop))
  allocate(Glob_NonlinParam(npt,Kstop))
  Glob_History(1:Glob_CurrBasisSize)=TempHistory(1:Glob_CurrBasisSize)
  Glob_FuncNum(1:Glob_CurrBasisSize)=TempFunc(1:Glob_CurrBasisSize)
  Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)=TempParam(1:npt,1:Glob_CurrBasisSize)
  deallocate(TempParam)
  deallocate(TempFunc)
  deallocate(TempHistory)
endif

!Allocate some global arrays
allocate(Glob_H(Kstop,Kstop))
allocate(Glob_S(Kstop,Kstop))
allocate(Glob_diagH(Kstop))
allocate(Glob_diagS(Kstop))
allocate(Glob_D(2*npt,Kstep,Kstop))
allocate(Glob_c(Kstop))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',Kstop,Kstop,Kstop,Kstop)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*Kstop 
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*Kstop))
allocate(Glob_IWorkForZHEGVX(5*Kstop))

!Allocate workspace for EnergyGB
allocate(Glob_WkGR(Kstep*npt))

!Allocate workspace
allocate(ParSet(npt,Kstep))
allocate(ParSetBest(npt,Kstep))
allocate(x(nvmax))
allocate(x_best(nvmax))
allocate(grad(nvmax))

!Allocate arrays used by DRMNG
nvmax=npt*Kstep
allocate(D(nvmax))
LV=71+nvmax*(nvmax+13)/2 + 1
allocate(V(LV))
allocate(V_init(LV))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
call DIVSET(ALG,IV_init,LIV,LV,V_init)
IV_init(17)=1000000 
IV_init(18)=1000000
IV_init(19)=0 !set summary print format
IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
V_init(31)=0.D0*ONE
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
IV_init(1)=12 !DIVSET has been called and some default values were changed

call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Calculating the initial energy
ErrCode=0
if (Kstart>1) then
  if (IsSwapFileOK) then
    if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
	Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.false.,ErrCode)
  else
    if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
	Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.true.,ErrCode)
  endif
else
    Glob_CurrEnergy=huge(Glob_CurrEnergy)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in BasisEnlG: initial energy cannot be computed'
  stop
endif

if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy


rgm1_counter=0
rgm2_counter=0
ms1=ZERO
ms2=ZERO
K=Kstart-1

!Main loop begins here
do while (K<Kstop)
  if (K+Kstep<=Kstop) then
	nfo=Kstep
	nfru=K    
	K=K+Kstep
  else
    nfo=Kstop-K
	nfru=K
	K=Kstop
  endif
  Glob_nfa=K
  Glob_nfru=nfru
  Glob_nfo=nfo
  nfrup1=nfru+1
  nv=nfo*npt
  E_init=Glob_CurrEnergy
  if (Glob_ProcID==0) then
    write(*,*)
    write(*,'(1x,a22,i5)') 'Current basis size is ',Glob_CurrBasisSize
	if (nfo>1) then
	  write(*,'(1x,a21,i5,a10,i5)') 'Selecting functions  ',nfrup1,'  through ',K
    else
	  write(*,'(1x,a19,i5)') 'Selecting function ',K
	endif
  endif
  
  IsOverlapBad=.true.
  IsAnyLinCoeffBad=.true.
  AttemptToGetGoodFunc=1
  do while ((IsOverlapBad.or.IsAnyLinCoeffBad).and. &
            (AttemptToGetGoodFunc<=Glob_BadOverlapOrLinCoeffLim))
    !Random selection 
    NumOfFailures=0
    IsEnergyImproved=.false.
	wbfu=0
	wmu=0
    do i=1,NTrials
      if (Glob_ProcID==0) call GenerateTrialParam(nfo,ParSet,wbfu_t,wmu_t)
	  call MPI_BCAST(ParSet,npt*nfo,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      Glob_NonlinParam(1:npt,nfrup1:K)=ParSet(1:npt,1:nfo)
	  Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
	  if (ErrCode==0) then
        if (Evalue<Glob_CurrEnergy) then
          Glob_CurrEnergy=Evalue
		  ParSetBest(1:npt,1:nfo)=ParSet(1:npt,1:nfo)
		  IsEnergyImproved=.true.
		  wbfu=wbfu_t
		  wmu=wmu_t
	    endif
	  else
        NumOfFailures=NumOfFailures+1
	  endif
    enddo
    if (NumOfFailures*ONE/NTrials>Glob_MaxFracOfTrialFailsAllowed) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlG: the limit of eigenvalue problem solution failures'
	    write(*,*) 'in random selection process is exceeded'
	    write(*,'(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
	                           (100*NumOfFailures*ONE)/NTrials,'%'
      endif
	  stop
    endif
    if (.not.(IsEnergyImproved)) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlG: random selection did not result'
	    write(*,*) 'in any energy improvement'
      endif
	  stop
    endif
    Glob_NonlinParam(1:npt,nfrup1:K)=ParSetBest(1:npt,1:nfo)
    Glob_CurrEnergy=EnergyGA(nfrup1,K,.true.,ErrCode)
    if (Glob_ProcID==0) then
	  write (*,'(1x,a2,d23.16,a23,i5)') 'E=',Glob_CurrEnergy,'  prototype function is',wbfu
      do i=1,nfo
        write(*,'(1x,i5,a2,<npt>d23.16)') nfru+i,': ',ParSetBest(1:npt,i)
      enddo
	  write (*,'(1x,a31)') 'Optimizing nonlinear parameters'		    
    endif
  
    !Optimization of the nonlinear parameters
    select case (OptimizationType)
    case(0) !No optimization 
	  do i=1,nfo
	    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  enddo    
    case(1)
    
	  !Setting IV and V values as was in their initial copies
	  IV(1:LIV)=IV_init(1:LIV)
	  V(1:LV)=V_init(1:LV)
 
	  nv=nfo*npt
	  do i=1,nfo
	    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  enddo
	  if (nfru>=nfo) then
	    t=max(abs((E_init-Glob_CurrEnergy))/(abs(E_init)+abs(Glob_CurrEnergy)),10000*epsilon(Glob_CurrEnergy))
	  else
	    t=ONE
	  endif 	  
      do i=1,nfo
        !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
        do j=1,np
          !Make sure none of the D(i) will be zero or smaller than the threshold
          !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
          D(npt*(i-1)+j)=t
        enddo
        !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
        D(npt*(i-1)+np+1:npt*i)=t
      enddo

	  ExitNeeded=.false.
	  NumOfFailures=0
      NumOfEnergyEval=0
      NumOfGradEval=0
      if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
      E_best=Glob_CurrEnergy
      x_best(1:nfo*npt)=x(1:nfo*npt)

      do while (.not.(ExitNeeded))  
	    if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
        call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
	    select case (IV(1))
        case (1) !Only energy is needed
          call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	      do i=1,nfo
	        Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      enddo
	      Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
          NumOfEnergyEval=NumOfEnergyEval+1
		  if (ErrCode/=0) then
            NumOfFailures=NumOfFailures+1
		    IV(2)=0
		  else
            Glob_CurrEnergy=Evalue
            if (Evalue<E_best) then 
              E_best=Evalue
              x_best(1:nfo*npt)=x(1:nfo*npt)
            endif  
		  endif
	    case (2) !Only gradient is needed
          call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	      do i=1,nfo
	        Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      enddo
          call EnergyGB(Evalue,grad,.true.,ErrCode)
          NumOfGradEval=NumOfGradEval+1
		  if (ErrCode/=0) then
            NumOfFailures=NumOfFailures+1
		    IV(2)=0
		  else
            if (Evalue<E_best) then 
              E_best=Evalue
              x_best(1:nfo*npt)=x(1:nfo*npt)
            endif 		      
		  endif
		  !===================================
		  !The lines above need be uncommented when finite
		  !difference gradient is shown 
		  !allocate(fx(nvmax))
          !allocate(fgrad(nvmax))	
		  !do j=1,nv
          !  fx(1:nv)=x(1:nv)
		  !  deltax=x(j)*3.0D-07
          !  fx(j)=x(j)+deltax
 	      !  do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
	      !  enddo
          !  Evalue1=EnergyGA(nfrup1,K,.true.,ErrCode)
          !  fx(j)=x(j)-deltax
 	      !  do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
	      !  enddo
          !  Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
		  !  fgrad(j)=(Evalue1-Evalue)/(2*deltax)
		  !enddo
		  !do j=1,nv
          !  write(*,*) j,'  ',grad(j),'  ' ,fgrad(j)
		  !enddo
          !write(*,*)
 	      !do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      !enddo
		  !deallocate(fgrad)
          !deallocate(fx)
		  !===================================	     
	    case (3:8) !Some kind of convergence has been reached
          ExitNeeded=.true.
        case (9:10)   !Function evaluation limit has been reached.
	                !This never suppose to happen because we
				    !count the number of function evaluations ourselves.
          ExitNeeded=.true.
	    endselect  
	    if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	      if (Glob_ProcID==0) then
            write(*,*) 'Error in BasisEnlG: number of failures in energy or gradient'
		    write(*,*) 'calculations during the optimization of nonlinear parameters'
		    write(*,*) 'exceeded limit' 
		  endif
	      stop
	    endif
	    if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
      enddo
    endselect !(OptimizationType)
    
    !Calculate the energy and linear coefficients at the best point found
    do i=1,nfo
	  Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
    enddo
!    Glob_CurrEnergy=EnergyGAM(nfrup1,K,.true.,ErrCode)
!
! call to EnergyGAM is replaced by a call to EnergyGA
!
    Glob_CurrEnergy=EnergyGA(nfrup1,K,.true.,ErrCode)
    !!We run EnergyGA again because EnergyGAM might give slightly different
    !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
    !Glob_CurrEnergy=EnergyGA(nfrup1,K,.false.,ErrCode)    
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlG: failed to evaluate energy after the optimization'
	    write(*,*) 'of nonlinear parameters'
	  endif 
	  stop     
    endif
	
	!Checking if overlaps are OK (only in case OverlapThreshold>ZERO)
    IsOverlapBad=.false.
    if (OverlapThreshold>ZERO) then
      ii=0
      do i=nfrup1,K
        do j=1,i-1
          if (abs(Glob_S(i,j))>OverlapThreshold) then
            ii=ii+1
			IsOverlapBad=.true.
	        if (Glob_ProcID==0) then
	          if (ii==1) then
                write(*,*) 'Warning: overlap of the following functions exceeds threshold'
		        write(*,*) 'Generated basis function(s) are rejected and a new attempt to'
		        write(*,*) 'generate them will be made'
		      endif
	          write(*,'(i5,a1,i5,i5,a7,d23.16,a1,d23.16,a1)') &
		           ii,':',i,j,'    S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),')'
            endif
	      endif
        enddo
      enddo
    endif
	!Checking if liniar coefficients are OK (only in case LinCoeffThreshold>ZERO)
	IsAnyLinCoeffBad=.false.
	if (LinCoeffThreshold>ZERO) then
      ii=0
	  do i=1,K
        if(abs(Glob_c(i))>LinCoeffThreshold) then
          ii=ii+1
		  IsAnyLinCoeffBad=.true.
		  if (Glob_ProcID==0) then
            if (ii==1) then
              write(*,*) 'Warning: absolute value of linear parameters of the'
			  write(*,*) 'following functions exceeds threshold'
			  write(*,*) 'Generated basis function(s) are rejected and a new attempt to'
		      write(*,*) 'generate them will be made'
			endif
			write(*,'(i5,a1,i5,a7,d23.16,a1,d23.16,a1)') &
			      ii,':',i,'    c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')'
		  endif
		endif
      enddo
	endif
	AttemptToGetGoodFunc=AttemptToGetGoodFunc+1

  enddo !(IsOverlapBad.or.IsAnyLinCoeffBad).and.(AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)
  
  if (Glob_ProcID==0) then
    write (*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
    write (*,'(1x,a2,d23.16)') 'E=',Glob_CurrEnergy
    do i=1,nfo
      write(*,'(1x,i5,a2,<npt>d23.16)') nfru+i,': ',Glob_NonlinParam(1:npt,nfru+i)
    enddo		    
  endif

  !getting statistics about the distribution of generated parameters 
  if (Glob_ProcID==0) then  
    if (wmu==1) then 
	  rgm1_counter=rgm1_counter+1
	  if (nfru>nfo) then
	    do i=1,nfo
          do j=1,npt
		    t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/Glob_NonlinParam(j,wbfu+i-1)
	        ms1=ms1+abs(t)
          enddo
	    enddo
	  endif
    endif
    if (wmu==2) then 
	  rgm2_counter=rgm2_counter+1
	  if (nfru>nfo) then
	    do i=1,nfo
          do j=1,npt
		    t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/Glob_NonlinParam(j,wbfu+i-1)
	        ms2=ms2+abs(t)        
		  enddo
	    enddo
	  endif
    endif
  endif

  Glob_CurrBasisSize=K
  do i=1,nfo
    Glob_History(nfru+i)%Energy=Glob_CurrEnergy
    Glob_History(nfru+i)%CyclesDone=0
	Glob_History(nfru+i)%InitFuncAtLastStep=0
    Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt=0
  enddo
  do i=1,nfo
    Glob_FuncNum(nfru+i)=nfru+i
  enddo

  if (Glob_ProcID==0) call SaveResults(Sort='no')
enddo
!Main loop ends here

call StoreMatricesInSwapFile()   

!Dellocate arrays used by DRMNG
deallocate(V_init)
deallocate(V)
deallocate(D)

!Deallocate workspace
deallocate(grad)
deallocate(x_best)
deallocate(x)
deallocate(ParSetBest)
deallocate(ParSet)

!Deallocate workspace for EnergyGB 
deallocate(Glob_WkGR)

!Deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX)
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!Deallocate global arrays 
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_D)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(GLob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) then
  write(*,*) 'Random selection statistics:'
  write(*,*) 'Method 1 of generating basis functions was used ',rgm1_counter,' times'
  if (rgm1_counter/=0) write(*,*) 'Average shift factor from prototype function is ',ms1/(npt*rgm1_counter)
  write(*,*) 'Method 2 of generating basis functions was used ',rgm2_counter,' times'
  if (rgm2_counter/=0) write(*,*) 'Average shift factor from prototype function is ',ms2/(npt*rgm2_counter)   
  write(*,*)
  write(*,'(1x,a30)') 'Routine BasisEnlG has finished'
endif

end subroutine BasisEnlG



subroutine BasisEnlI(Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval, &
OverlapThreshold,LinCoeffThreshold) 
!Subroutine BasisEnlI enlarges the basis set from initial
!Kstart-1 functions to Kstop functions by means of trials of 
!randomly selected candidates and then optimizing them. At each 
!step it generates a set of Kstep functions for NTrials times
!based on the existing distribution of nonlinear papameters. 
!Only the set that lowers the energy the most is left. Then this
!set is optimized according to the value of OptimizationType.
!MaxEnergyEval is the limit of the energy evaluations for this
!optimization process. 
!Each newly accepted basis function is checked for pair linear
!dependency with other functions in the basis. If the absolute
!value of the overlap is greater than OverlapThreshold then
!such function is rejected and random trial and optimization 
!take place again until either a good function/functions are 
!generated or a certain limit of such repetitions is reached 
!(the limit is given by the value of global constant 
!Glob_BadOverlapOrLinCoeffLim). If the value of OverlapThreshold 
!passed is negative or zero then no such check is performed. 
!A Similar check is performed for all linear parameters. If
!adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions 
!are rejected and the procedure is repeted until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avaid this check it is enough to set LinCoeffThreshold to a value 
!equal or smaller than zero.
!It is assumed that the basis of Kstart-1 functions has already
!been selected and the corresponding nonlinear parameters are
!stored in proper global array. Subroutines EnergyIA, EnergyIB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)     :: Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval
real(dprec),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
integer      i,j,K,AttemptToGetGoodFunc,ii
integer      np,npt,nfo,nfa,nfru,nfrup1,nvmax,nv
integer      OpenFileErr,ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
logical      IsSwapFileOK,IsEnergyImproved,ExitNeeded
logical      IsOverlapBad,IsAnyLinCoeffBad
integer      wbfu_t,wmu_t,wbfu,wmu,rgm1_counter,rgm2_counter
real(dprec)  ms1,ms2
real(dprec)  Evalue,E_init,E_best
real(dprec)  t
real(dprec),allocatable,dimension(:,:)   :: ParSet,ParSetBest
real(dprec),allocatable,dimension(:)     :: x,x_best,grad 
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V,V_init
integer,parameter    :: LIV=60
integer                 IV(LIV), IV_init(LIV)
integer                 LV
integer                 ALG
!Allocatable work space
real(dprec),allocatable,dimension(:)            :: WorkBuffReal
integer,allocatable,dimension(:)                :: WorkBuffInt
type(Glob_HistoryStep),allocatable,dimension(:) :: TempHistory
real(dprec),allocatable,dimension(:,:)          :: TempParam
integer,allocatable,dimension(:)                :: TempFunc
!==================================================== 
!These variables are used when a finite difference gradient is computed               
!real(dprec),allocatable,dimension(:)     ::    fx,fgrad
!real(dprec)                                    deltax,Evalue1
!====================================================

if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine BasisEnlI has started'
  write(*,'(1x,a7,i5)') 'Kstart=',Kstart
  write(*,'(1x,a7,i5)') 'Kstop= ',Kstop
  write(*,'(1x,a7,i5)') 'Kstep= ',Kstep
  write(*,'(1x,a17,i4)') 'OptimizationType=',OptimizationType
  write(*,'(1x,a14,i7)') 'MaxEnergyEval=',MaxEnergyEval
endif

!Setting the values of some global variables
Glob_GSEPSolutionMethod='I'
Glob_OverlapPenaltyAllowed=.false.
Glob_nfa=Kstart+Kstep
Glob_nfo=Kstep
Glob_HSLeadDim=Kstop
Glob_HSBuffLen=Kstop*Kstep
np=Glob_np
npt=Glob_npt
nfo=Glob_nfo
nfa=Glob_nfa
nvmax=Kstep*Glob_npt

!Reallocate arrays that contain the information
!about basis functions and optimization process.

if (Glob_UseReallocFile) then
  !Temporarily store the information in a file
  if (Glob_ProcID==0) then
    if (Glob_CurrBasisSize>0) then
      write(*,'(1x,a47)',advance='no') 'Reallocating some arrays using external file...'
      open(1,file=Glob_ReallocFileName,form='binary',status='replace')
      write(1) Glob_History(1:Glob_CurrBasisSize)
      write(1) Glob_FuncNum(1:Glob_CurrBasisSize)
      write(1) Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
      close(1)
	endif
  endif
  deallocate(Glob_NonlinParam)
  deallocate(Glob_FuncNum)
  deallocate(Glob_History)
  allocate(Glob_History(Kstop))
  allocate(Glob_FuncNum(Kstop))
  allocate(Glob_NonlinParam(npt,Kstop))
  if (Glob_ProcID==0) then
    if (Glob_CurrBasisSize>0) then
      open(1,file=Glob_ReallocFileName,form='binary',status='old',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        read(1) Glob_History(1:Glob_CurrBasisSize)
        read(1) Glob_FuncNum(1:Glob_CurrBasisSize)
        read(1,iostat=OpenFileErr) Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
      endif
      close(1) 
      call MPI_BCAST(OpenFileErr,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  endif
  if ((OpenFileErr/=0).and.(Glob_CurrBasisSize>0)) then
    if (Glob_ProcID==0) then
	  write(*,*)
	  write(*,*) 'Error in BasisEnlI: cannot read data from ',Glob_ReallocFileName
    endif
    stop
  endif
  if (Glob_ProcID==0) then
    open(1,file=Glob_ReallocFileName,form='binary',status='replace',iostat=OpenFileErr)
    write(1)  'This temporary file is empty'
    close(1)
  endif
  allocate(WorkBuffReal(Glob_CurrBasisSize))
  allocate(WorkBuffInt(Glob_CurrBasisSize))
  do i=1,Glob_CurrBasisSize
    WorkBuffReal(i)=Glob_History(i)%Energy
  enddo
  call MPI_BCAST(WorkBuffReal,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%Energy=WorkBuffReal(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%CyclesDone
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%CyclesDone=WorkBuffInt(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%InitFuncAtLastStep
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%InitFuncAtLastStep=WorkBuffInt(i)
  enddo
  do i=1,Glob_CurrBasisSize
    WorkBuffInt(i)=Glob_History(i)%NumOfEnergyEvalDuringFullOpt
  enddo
  call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do i=1,Glob_CurrBasisSize
    Glob_History(i)%NumOfEnergyEvalDuringFullOpt=WorkBuffInt(i)
  enddo
  deallocate(WorkBuffReal)
  deallocate(WorkBuffInt)
  call MPI_BCAST(Glob_NonlinParam,npt*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_ProcID==0) then
    write(*,*) 'completed'
  endif
else ! if (Glob_UseReallocFile)
  allocate(TempHistory(Glob_CurrBasisSize))
  allocate(TempFunc(Glob_CurrBasisSize))
  allocate(TempParam(npt,Glob_CurrBasisSize))
  TempHistory(1:Glob_CurrBasisSize)=Glob_History(1:Glob_CurrBasisSize)
  TempFunc(1:Glob_CurrBasisSize)=Glob_FuncNum(1:Glob_CurrBasisSize)
  TempParam(1:npt,1:Glob_CurrBasisSize)=Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)
  deallocate(Glob_NonlinParam)
  deallocate(Glob_FuncNum)
  deallocate(Glob_History)
  allocate(Glob_History(Kstop))
  allocate(Glob_FuncNum(Kstop))
  allocate(Glob_NonlinParam(npt,Kstop))
  Glob_History(1:Glob_CurrBasisSize)=TempHistory(1:Glob_CurrBasisSize)
  Glob_FuncNum(1:Glob_CurrBasisSize)=TempFunc(1:Glob_CurrBasisSize)
  Glob_NonlinParam(1:npt,1:Glob_CurrBasisSize)=TempParam(1:npt,1:Glob_CurrBasisSize)
  deallocate(TempParam)
  deallocate(TempFunc)
  deallocate(TempHistory)
endif

!Allocate some global arrays
allocate(Glob_H(Kstop,Kstop))
allocate(Glob_S(Kstop,Kstop))
allocate(Glob_diagS(Kstop))
allocate(Glob_invD(Kstop))
allocate(Glob_D(2*npt,Kstep,Kstop))
allocate(Glob_c(Kstop))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GHEPIIS, which is called
!inside EnergyIA, EnergyIAM, and EnergyIB
allocate(Glob_WorkForGHEPIIS(Kstop))
allocate(Glob_LastEigvector(Kstop))
Glob_LastEigvector(1:Kstop)=ONE

!Allocate workspace for EnergyIB
allocate(Glob_WkGR(Kstep*npt))

!Allocate workspace
allocate(ParSet(npt,Kstep))
allocate(ParSetBest(npt,Kstep))
allocate(x(nvmax))
allocate(x_best(nvmax))
allocate(grad(nvmax))

!Allocate arrays used by DRMNG
nvmax=npt*Kstep
allocate(D(nvmax))
LV=71+nvmax*(nvmax+13)/2 + 1
allocate(V(LV))
allocate(V_init(LV))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
call DIVSET(ALG,IV_init,LIV,LV,V_init)
IV_init(17)=1000000 
IV_init(18)=1000000
IV_init(19)=0 !set summary print format
IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
V_init(31)=0.D0*ONE
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
IV_init(1)=12 !DIVSET has been called and some default values were changed

call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Calculating the initial energy
call linalg_setparam(Glob_CurrBasisSize)
ErrCode=0
if (Kstart>1) then
  if (IsSwapFileOK) then
    if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
	Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.false.,ErrCode)
  else
    if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
	Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.true.,ErrCode)
  endif
else
    Glob_CurrEnergy=huge(Glob_CurrEnergy)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in BasisEnlI: initial energy cannot be computed'
  stop
endif

if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy


rgm1_counter=0
rgm2_counter=0
ms1=ZERO
ms2=ZERO
K=Kstart-1

!Main loop begins here
do while (K<Kstop)
  if (K+Kstep<=Kstop) then
	nfo=Kstep
	nfru=K    
	K=K+Kstep
  else
    nfo=Kstop-K
	nfru=K
	K=Kstop
  endif
  Glob_nfa=K
  Glob_nfru=nfru
  Glob_nfo=nfo
  nfrup1=nfru+1
  nv=nfo*npt
  E_init=Glob_CurrEnergy
  Glob_InvItTempCounter1=0
  Glob_InvItTempCounter2=0  
  if (Glob_ProcID==0) then
    write(*,*)
    write(*,'(1x,a22,i5)') 'Current basis size is ',Glob_CurrBasisSize
	if (nfo>1) then
	  write(*,'(1x,a21,i5,a10,i5)') 'Selecting functions  ',nfrup1,'  through ',K
    else
	  write(*,'(1x,a19,i5)') 'Selecting function ',K
	endif
  endif
  call linalg_setparam(K)
  IsOverlapBad=.true.
  IsAnyLinCoeffBad=.true.
  AttemptToGetGoodFunc=1
  do while ((IsOverlapBad.or.IsAnyLinCoeffBad).and. &
            (AttemptToGetGoodFunc<=Glob_BadOverlapOrLinCoeffLim))
    !Random selection 
    NumOfFailures=0
    IsEnergyImproved=.false.
	wbfu=0
	wmu=0
    do i=1,NTrials
      if (Glob_ProcID==0) call GenerateTrialParam(nfo,ParSet,wbfu_t,wmu_t)
	  call MPI_BCAST(ParSet,npt*nfo,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      Glob_NonlinParam(1:npt,nfrup1:K)=ParSet(1:npt,1:nfo)
	  Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
	  if (ErrCode==0) then
        if (Evalue<Glob_CurrEnergy) then
          Glob_CurrEnergy=Evalue
		  ParSetBest(1:npt,1:nfo)=ParSet(1:npt,1:nfo)
		  IsEnergyImproved=.true.
		  wbfu=wbfu_t
		  wmu=wmu_t
	    endif
	  else
        NumOfFailures=NumOfFailures+1
	  endif
    enddo
    if (NumOfFailures*ONE/NTrials>Glob_MaxFracOfTrialFailsAllowed) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlI: the limit of eigenvalue problem solution failures'
	    write(*,*) 'in random selection process is exceeded'
	    write(*,'(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
	                           (100*NumOfFailures*ONE)/NTrials,'%'
      endif
	  stop
    endif
    if (.not.(IsEnergyImproved)) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlI: random selection did not result'
	    write(*,*) 'in any energy improvement'
      endif
	  stop
    endif
    Glob_NonlinParam(1:npt,nfrup1:K)=ParSetBest(1:npt,1:nfo)
    Glob_CurrEnergy=EnergyIA(nfrup1,K,.true.,ErrCode)
    if (Glob_ProcID==0) then
	  write (*,'(1x,a2,d23.16,a23,i5)') 'E=',Glob_CurrEnergy,'  prototype function is',wbfu
      do i=1,nfo
        write(*,'(1x,i5,a2,<npt>d23.16)') nfru+i,': ',ParSetBest(1:npt,i)
      enddo
	  write (*,'(1x,a31)') 'Optimizing nonlinear parameters'		    
    endif
  
    !Optimization of the nonlinear parameters
    select case (OptimizationType)
    case(0) !No optimization 
	  do i=1,nfo
	    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  enddo    
    case(1)
    
	  !Setting IV and V values as was in their initial copies
	  IV(1:LIV)=IV_init(1:LIV)
	  V(1:LV)=V_init(1:LV)
 
	  nv=nfo*npt
	  do i=1,nfo
	    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  enddo
	  if (nfru>=nfo) then
	    t=max(abs((E_init-Glob_CurrEnergy))/(abs(E_init)+abs(Glob_CurrEnergy)),10000*epsilon(Glob_CurrEnergy))
	  else
	    t=ONE
	  endif 	  
      do i=1,nfo
        !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
        do j=1,np
          !Make sure none of the D(i) will be zero or smaller than the threshold
          !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
          D(npt*(i-1)+j)=t
        enddo
        !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
        D(npt*(i-1)+np+1:npt*i)=t
      enddo

	  ExitNeeded=.false.
	  NumOfFailures=0
      NumOfEnergyEval=0
      NumOfGradEval=0
      if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
      E_best=Glob_CurrEnergy
      x_best(1:nfo*npt)=x(1:nfo*npt)

      do while (.not.(ExitNeeded))  
	    if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
        call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
	    select case (IV(1))
        case (1) !Only energy is needed
          call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	      do i=1,nfo
	        Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      enddo
	      Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
          NumOfEnergyEval=NumOfEnergyEval+1
		  if (ErrCode/=0) then
            NumOfFailures=NumOfFailures+1
		    IV(2)=0
		  else
            Glob_CurrEnergy=Evalue
            if (Evalue<E_best) then 
              E_best=Evalue
              x_best(1:nfo*npt)=x(1:nfo*npt)
            endif  
		  endif
	    case (2) !Only gradient is needed
          call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	      do i=1,nfo
	        Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      enddo
          call EnergyIB(Evalue,grad,.true.,ErrCode)
          NumOfGradEval=NumOfGradEval+1
		  if (ErrCode/=0) then
            NumOfFailures=NumOfFailures+1
		    IV(2)=0
		  else
            if (Evalue<E_best) then 
              E_best=Evalue
              x_best(1:nfo*npt)=x(1:nfo*npt)
            endif 		      
		  endif
		  !===================================
		  !The lines above need be uncommented when finite
		  !difference gradient is shown 
		  !allocate(fx(nvmax))
          !allocate(fgrad(nvmax))	
		  !do j=1,nv
          !  fx(1:nv)=x(1:nv)
		  !  deltax=x(j)*3.0D-07
          !  fx(j)=x(j)+deltax
 	      !  do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
	      !  enddo
          !  Evalue1=EnergyGA(nfrup1,K,.true.,ErrCode)
          !  fx(j)=x(j)-deltax
 	      !  do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
	      !  enddo
          !  Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
		  !  fgrad(j)=(Evalue1-Evalue)/(2*deltax)
		  !enddo
		  !do j=1,nv
          !  write(*,*) j,'  ',grad(j),'  ' ,fgrad(j)
		  !enddo
          !write(*,*)
 	      !do i=1,nfo
	      !    Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	      !enddo
		  !deallocate(fgrad)
          !deallocate(fx)
		  !===================================	     
	    case (3:8) !Some kind of convergence has been reached
          ExitNeeded=.true.
        case (9:10)   !Function evaluation limit has been reached.
	                !This never suppose to happen because we
				    !count the number of function evaluations ourselves.
          ExitNeeded=.true.
	    endselect  
	    if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	      if (Glob_ProcID==0) then
            write(*,*) 'Error in BasisEnlI: number of failures in energy or gradient'
		    write(*,*) 'calculations during the optimization of nonlinear parameters'
		    write(*,*) 'exceeded limit' 
		  endif
	      stop
	    endif
	    if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
      enddo
    endselect !(OptimizationType)
    
    !Calculate the energy and linear coefficients at the best point found
    do i=1,nfo
	  Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
    enddo
    Glob_CurrEnergy=EnergyIAM(nfrup1,K,.true.,ErrCode)    
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Error in BasisEnlI: failed to evaluate energy after the optimization'
	    write(*,*) 'of nonlinear parameters'
	  endif 
	  stop     
    endif
	
	!Checking if overlaps are OK (only in case OverlapThreshold>ZERO)
    IsOverlapBad=.false.
    if (OverlapThreshold>ZERO) then
      ii=0
      do i=nfrup1,K
        do j=1,i-1
          if (abs(Glob_S(i,j))>OverlapThreshold) then
            ii=ii+1
			IsOverlapBad=.true.
	        if (Glob_ProcID==0) then
	          if (ii==1) then
                write(*,*) 'Warning: overlap of the following functions exceeds threshold'
		        write(*,*) 'Generated basis function(s) are rejected and a new attempt to'
		        write(*,*) 'generate them will be made'
		      endif
	          write(*,'(i5,a1,i5,i5,a7,d23.16,a1,d23.16,a1)') &
		           ii,':',i,j,'    S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),')'
            endif
	      endif
        enddo
      enddo
    endif
	!Checking if liniar coefficients are OK (only in case LinCoeffThreshold>ZERO)
	IsAnyLinCoeffBad=.false.
	if (LinCoeffThreshold>ZERO) then
      ii=0
	  do i=1,K
        if(abs(Glob_c(i))>LinCoeffThreshold) then
          ii=ii+1
		  IsAnyLinCoeffBad=.true.
		  if (Glob_ProcID==0) then
            if (ii==1) then
              write(*,*) 'Warning: absolute value of linear parameters of the'
			  write(*,*) 'following functions exceeds threshold'
			  write(*,*) 'Generated basis function(s) are rejected and a new attempt to'
		      write(*,*) 'generate them will be made'
			endif
			write(*,'(i5,a1,i5,a7,d23.16,a1,d23.16,a1)') &
			      ii,':',i,'    c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')'
		  endif
		endif
      enddo
	endif
	AttemptToGetGoodFunc=AttemptToGetGoodFunc+1

  enddo !(IsOverlapBad.or.IsAnyLinCoeffBad).and.(AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)
  
  if (Glob_ProcID==0) then
    write (*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
    write (*,'(1x,a2,d23.16)') 'E=',Glob_CurrEnergy
    do i=1,nfo
      write(*,'(1x,i5,a2,<npt>d23.16)') nfru+i,': ',Glob_NonlinParam(1:npt,nfru+i)
    enddo		    
  endif

  !getting statistics about the distribution of generated parameters 
  if (Glob_ProcID==0) then
    if (wmu==1) then 
	  rgm1_counter=rgm1_counter+1
	  if (nfru>nfo) then
	    do i=1,nfo
          do j=1,npt
		    t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/Glob_NonlinParam(j,wbfu+i-1)
	        ms1=ms1+abs(t)
          enddo
	    enddo
	  endif
    endif
    if (wmu==2) then 
	  rgm2_counter=rgm2_counter+1
	  if (nfru>nfo) then
	    do i=1,nfo
          do j=1,npt
		    t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/Glob_NonlinParam(j,wbfu+i-1)
	        ms2=ms2+abs(t)        
		  enddo
	    enddo
	  endif
    endif
  endif
 
  Glob_CurrBasisSize=K
  do i=1,nfo
    Glob_History(nfru+i)%Energy=Glob_CurrEnergy
    Glob_History(nfru+i)%CyclesDone=0
	Glob_History(nfru+i)%InitFuncAtLastStep=0
    Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt=0
  enddo
  do i=1,nfo
    Glob_FuncNum(nfru+i)=nfru+i
  enddo

  if (Glob_ProcID==0) call SaveResults(Sort='no')
enddo
!Main loop ends here

call StoreMatricesInSwapFile()   

!Dellocate arrays used by DRMNG
deallocate(V_init)
deallocate(V)
deallocate(D)

!Deallocate workspace
deallocate(grad)
deallocate(x_best)
deallocate(x)
deallocate(ParSetBest)
deallocate(ParSet)

!Deallocate workspace for EnergyGB 
deallocate(Glob_WkGR)

!Deallocate workspace for subroutine GHEPIIS
deallocate(Glob_LastEigvector)
deallocate(Glob_WorkForGHEPIIS)

!Deallocate global arrays 
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_D)
deallocate(Glob_invD)
deallocate(Glob_diagS)
deallocate(GLob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) then
  write(*,*) 'Random selection statistics:'
  write(*,*) 'Method 1 of generating basis functions was used ',rgm1_counter,' times'
  if (rgm1_counter/=0) write(*,*) 'Average shift factor from prototype function is ',ms1/(npt*rgm1_counter)
  write(*,*) 'Method 2 of generating basis functions was used ',rgm2_counter,' times'
  if (rgm2_counter/=0) write(*,*) 'Average shift factor from prototype function is ',ms2/(npt*rgm2_counter)   
  write(*,*)
  write(*,'(1x,a30)') 'Routine BasisEnlI has finished'
endif

end subroutine BasisEnlI



subroutine OptCycleG(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                     NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold) 
!Subroutine OptCycleG performs an optimization of the basis set of K functions 
!by means of cyclic optimization of a set of NumOfFuncToOpt functions at a time.
!After a step is made, the subroutine takes another set of NumOfFuncToOpt 
!functions, whose numbers are shifted by NumFuncToShift (normally, one
!would probably want to use NumFuncToShift=NumOfFuncToOpt). This strategy
!is applied NumCycles times in a row to functions whose numbers range from 
!FuncBegin to FuncEnd. MaxEnergyEval is the limit of the energy evaluations 
!for each step. If it happens that at some step the basis functions being
!optimized become linerly dependent with other basis functions in the basis
!then such a change is rejected and this step is omitted. Pair linear 
!dependency is checked by comparing the absolute value of the overlap
!with OverlapThreshold. If OverlapThreshold<0 then no linear dependency
!check is performed. A similar check is performed for all linear parameters. 
!If adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions 
!are rejected and the procedure is repeted until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avaid this check it is enough to set LinCoeffThreshold to a value 
!equal or smaller than zero. Subroutines EnergyGA, EnergyGB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)     :: K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift
integer,intent(in)     :: NumCycles, MaxEnergyEval
real(dprec),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
integer      i,j,m,ip,AttemptToGetGoodOverlap,ii
integer      np,npt,nfo,nfa,nfru,nfrup1,nv,nvmax,nfco,fbn,cbs,nfs
integer      CurrCycle,CurrFunc,CurrFuncBegin
integer      q,nr,OptIterCounter
integer      ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
integer      BlockSizeForZHEGVX
real(dprec)  Evalue,E_best
logical      IsSwapFileOK,ExitNeeded,LastIter
logical      IsOverlapBad,IsAnyLinCoeffBad
real(dprec)  t
real(dprec),allocatable,dimension(:,:) :: NonlinParamTemp
integer,allocatable,dimension(:)       :: FuncNumTemp
real(dprec),allocatable,dimension(:)   :: TempR
complex(dprec),allocatable,dimension(:):: TempC
real(dprec),allocatable,dimension(:)   :: x,grad,x_init,x_best
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V,V_init
integer,parameter    :: LIV=60
integer                 IV(LIV),IV_init(LIV)
integer                 LV
integer                 ALG

!Setting the values of some global variables
cbs=Glob_CurrBasisSize
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
np=Glob_np
npt=Glob_npt
nvmax=NumOfFuncToOpt*npt
Glob_HSLeadDim=Glob_CurrBasisSize
Glob_HSBuffLen=Glob_CurrBasisSize*NumOfFuncToOpt
Glob_nfa=Glob_CurrBasisSize
nfco=FuncEnd-FuncBegin+1
fbn=FuncBegin+Glob_CurrBasisSize-FuncEnd

!Checking if cyclic optimization is already completed for this basis size
if ((Glob_History(cbs)%CyclesDone>=NumCycles).and. &
    (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd)) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine OptCycleG has started'
  write(*,'(1x,a14,i5)') 'Basis size is ',cbs
  write(*,'(1x,a38)') 'Cyclic optimization of basis functions'
  write(*,'(1x,i5,a9,i5,a21)') FuncBegin,'  through',FuncEnd,' is already completed'
  write(*,'(1x,a20)') 'Exiting OptCycleG...'
  write (*,'(1x,a30)') 'Routine OptCycleG has finished'
  return
endif

if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine OptCycleG has started'
  write(*,'(1x,a14,i5)') 'Basis size is ',cbs
  write(*,'(1x,a38)') 'Cyclic optimization of basis functions'
  write(*,'(1x,i5,a9,i5,a19)') FuncBegin,'  through',FuncEnd,'  will be performed'
  write(*,'(1x,a14,i7)') 'MaxEnergyEval=',MaxEnergyEval
endif

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_D(2*npt,NumOfFuncToOpt,cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs)) 

!Allocate workspace for EnergyGB
allocate(Glob_WkGR(NumOfFuncToOpt*npt))

!Allocate workspace for SaveResults (it uses Sort='yes'
!option, which requires workspace)
allocate(Glob_IntWorkArrForSaveResults(cbs))

!Allocate arrays used by DRMNG
allocate(D(nvmax))
LV=71+nvmax*(nvmax+13)/2 + 1
allocate(V(LV))
allocate(V_init(LV))

!Allocate workspace
allocate(x(nvmax))
allocate(x_init(nvmax))
allocate(x_best(nvmax))
allocate(grad(nvmax))
allocate(NonlinParamTemp(npt,cbs-FuncBegin+1))
allocate(FuncNumTemp(cbs-FuncBegin+1))
allocate(TempR(cbs-FuncBegin+1))
allocate(TempC(cbs-FuncBegin+1))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again
	
!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
call DIVSET(ALG,IV_init,LIV,LV,V_init)
IV_init(17)=1000000 
IV_init(18)=1000000
IV_init(19)=0 !set summary print format
IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
V_init(31)=0.D0*ONE
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
IV_init(1)=12 !DIVSET has been called and some default values were changed

!Read swap file (if necessary) and distribute the data
call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Changing the order of the functions to be optimized and, if necessary, 
!the corresponding matrix elements to reverse
call ReverseFuncOrder(FuncBegin,FuncEnd)
if (IsSwapFileOK) call ReverseMatElemOrder(FuncBegin,FuncEnd)

!Shifting the set of basis functions to be optimized to
!the very end and, if necessary, doing the permutation of matrix elements
!to reflect this change.
call PermuteFunctions(FuncBegin,FuncEnd,FuncNumTemp,NonlinParamTemp)
if (IsSwapFileOK) call PermuteMatrixElements(FuncBegin,FuncEnd,TempR,TempC)

!If the last optimized function is greater than FuncEnd-1 or smaller than FuncBegin
!we need to change it so that the optimization begins from FuncBegin
if (Glob_History(cbs)%InitFuncAtLastStep<FuncBegin) Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
if (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd) then
  Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
  Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
endif

!We need to make an initial permutation to shift the functions that were 
!already optimized (if any) in the current optimization cycle
ip=Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
if (ip>0) then
  call PermuteFunctions(fbn,cbs-ip,FuncNumTemp,NonlinParamTemp)
  if (IsSwapFileOK) call PermuteMatrixElements(fbn,cbs-ip,TempR,TempC)   
endif

!Calculating the initial energy
if (IsSwapFileOK) then
  !Getting initial energy
  if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyGA(1,cbs,.false.,ErrCode)
else
  !Getting initial energy
  if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyGA(1,cbs,.true.,ErrCode)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in OptCycleG: initial energy cannot be computed'
  stop
endif

if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy


!Here comes main optimization cycle
do CurrCycle=Glob_History(cbs)%CyclesDone+1,NumCycles
  !Doing cycle number CurrCycle
  if (Glob_ProcID==0) then
    write(*,*)  
    write(*,'(1x,a5,i4,a10)') 'Cycle',CurrCycle,' has begun'
  endif
  CurrFuncBegin=Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
  q=FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
  OptIterCounter=0
  do CurrFunc=CurrFuncBegin,FuncEnd,NumOfFuncToShift
    !Note that CurrFunc counts functions as they were not shifted back
	!by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
	!initial input file.
    nfo=min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
	Glob_nfo=nfo	
    nv=nfo*npt
	nfru=cbs-nfo
	Glob_nfru=nfru
	nfrup1=nfru+1
	OptIterCounter=OptIterCounter+1
    if (Glob_ProcID==0) then
      write(*,*)
	  if (nfo>1) then
	    write(*,'(1x,a21,i5,a10,i5)') 'Optimizing functions ',CurrFunc,'  through ',CurrFunc+nfo-1
      else
	    write(*,'(1x,a21,i5)') 'Optimizing function ',CurrFunc
	  endif
    endif   
	
    !Permute basis functions so that the last nfo functions, which have just
	!been optimized are replaced with those that will be optimized at next step	 
	nfs=min(NumOfFuncToShift,nfo)
    if (OptIterCounter/=1) then
      i=cbs-nfo-nfs+1
	  j=cbs-nfs
      call PermuteFunctions(i,j,FuncNumTemp,NonlinParamTemp)
	  call PermuteMatrixElements(i,j,TempR,TempC)
	  nr=NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)
	  if (q-nfo>=nr*2) then
        !normal permutation, there is sufficient number of functions left
		m=cbs-nfo
		j=m-nr
		i=j-nr+1       
      else
        !abnormal permutation, not sufficient number of functions left
		m=cbs-nfo
		j=m-nr
		i=cbs-q+1
	  endif
      call PermuteFunctions2(i,j,m,FuncNumTemp,NonlinParamTemp)
	  call PermuteMatrixElements2(i,j,m,TempR,TempC) 
    endif

    if (Glob_ProcID==0) then
	  if (Glob_AreParamPrintedInCycleOptX) then
        write (*,*) 'Nonlinear parameters before optimization:'
        do i=1,nfo
          write(*,'(1x,i5,a2,<npt>d23.16)') Glob_FuncNum(nfru+i),': ',Glob_NonlinParam(1:npt,nfru+i)
        enddo	
	  endif		    
    endif

    !Setting IV and V values as was in their initial copies
	IV(1:LIV)=IV_init(1:LIV)
	V(1:LV)=V_init(1:LV)
	do i=1,nfo
	  x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  x_init((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	enddo
	
	t=max(ONE/(cbs*cbs*sqrt(ONE*cbs)),10000*epsilon(Glob_CurrEnergy))
    do i=1,nfo
      !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
      do j=1,np
        !Make sure none of the D(i) will be zero or smaller than the threshold
        !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        D(npt*(i-1)+j)=t
      enddo
      !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
      D(npt*(i-1)+np+1:npt*i)=t
    enddo

	ExitNeeded=.false.
	NumOfFailures=0
    NumOfEnergyEval=0
    NumOfGradEval=0
    if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
    E_best=Glob_CurrEnergy
    x_best(1:nfo*npt)=x(1:nfo*npt)    

    do while (.not.(ExitNeeded))  
	  if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
      call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
	  select case (IV(1))
      case (1) !Only energy is needed
        call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	    do i=1,nfo
	      Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	    enddo
	    Evalue=EnergyGA(nfrup1,cbs,.true.,ErrCode)
        NumOfEnergyEval=NumOfEnergyEval+1
		if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
		  IV(2)=0
		else
          Glob_CurrEnergy=Evalue
          if (Evalue<E_best) then 
            E_best=Evalue
            x_best(1:nfo*npt)=x(1:nfo*npt)
          endif           
		endif
	  case (2) !Only gradient is needed
        call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	    do i=1,nfo
	      Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	    enddo
        call EnergyGB(Evalue,grad,.true.,ErrCode)
        NumOfGradEval=NumOfGradEval+1
		if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
		  IV(2)=0
		else
          if (Evalue<E_best) then 
            E_best=Evalue
            x_best(1:nfo*npt)=x(1:nfo*npt)
          endif 			  
		endif	     
	  case (3:8) !Some kind of convergence has been reached
        ExitNeeded=.true.
      case (9:10) !Function evaluation limit has been reached.
	              !This is never supposed to happen because we
				  !count the number of function evaluations ourselves.
        ExitNeeded=.true.
	  endselect  
	  if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	    if (Glob_ProcID==0) then
          write(*,*) 'Error in OptCycleG: number of failures in energy or gradient'
		  write(*,*) 'calculations during the optimization of nonlinear parameters'
		  write(*,*) 'exceeded limit' 
		endif
	    stop
	  endif
	  if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
    enddo !while

    !Compute the energy and the linear coefficients at the best point found
    do i=1,nfo
	  Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
    enddo
!    Glob_CurrEnergy=EnergyGAM(nfrup1,cbs,.true.,ErrCode)
!
! call to EnergyGAM is replaced by a call to EnergyGA
!
    Glob_CurrEnergy=EnergyGA(nfrup1,cbs,.true.,ErrCode)
    !!We run EnergyGA again because EnergyGAM might give slightly different
    !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
    !Glob_CurrEnergy=EnergyGA(nfrup1,cbs,.false.,ErrCode)    
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Warning in OptCycleG: failed to evaluate energy after optimization'
	    write(*,*) 'of nonlinear parameters. The values of the nonlinear parameters'
		write(*,*) 'will be left unchanged.'
	  endif 
    endif
 
 	!Checking if overlap is OK (only in case OverlapThreshold>ZERO)
    IsOverlapBad=.false.
    if (OverlapThreshold>ZERO) then
      ii=0
      do i=nfrup1,cbs
        do j=1,i-1
          if (abs(Glob_S(i,j))>OverlapThreshold) then
            ii=ii+1
			IsOverlapBad=.true.
	        if (Glob_ProcID==0) then
	          if (ii==1) then
                write(*,*) 'Warning: overlap of the following functions exceeds threshold'
		        write(*,*) 'Nonlinear parameters will be left unchanged'
		      endif
	          write(*,'(i5,a1,i5,i5,a7,d23.16,a1,d23.16,a1)') &
		           ii,':',i,j,'    S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),')'
            endif
	      endif
        enddo
      enddo
    endif
	!Checking if liniar coefficients are OK (only in case LinCoeffThreshold>ZERO)
	IsAnyLinCoeffBad=.false.
	if (LinCoeffThreshold>ZERO) then
      ii=0
	  do i=1,cbs
        if(abs(Glob_c(i))>LinCoeffThreshold) then
          ii=ii+1
		  IsAnyLinCoeffBad=.true.
		  if (Glob_ProcID==0) then
            if (ii==1) then
              write(*,*) 'Warning: absolute value of linear parameters of the'
			  write(*,*) 'following functions exceeds threshold'
		      write(*,*) 'Nonlinear parameters will be left unchanged'
			endif
			write(*,'(i5,a1,i5,a7,d23.16,a1,d23.16,a1)') &
			      ii,':',i,'    c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')'
		  endif
		endif
      enddo
	endif

    if ((Glob_ProcID==0).and.(ErrCode==0).and.(.not.IsOverlapBad).and.(.not.IsAnyLinCoeffBad)) then
      write (*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
      write (*,'(1x,a2,d23.16)') 'E=',Glob_CurrEnergy
      if (Glob_AreParamPrintedInCycleOptX) then
        write (*,*) 'Nonlinear parameters after optimization:'
        do i=1,nfo
          write(*,'(1x,i5,a2,<npt>d23.16)') Glob_FuncNum(nfru+i),': ',Glob_NonlinParam(1:npt,nfru+i)
        enddo		    
      endif
    endif

    if ((ErrCode/=0).or.IsOverlapBad.or.IsAnyLinCoeffBad) then
	  !restore initial values of nonlinear parameters
      do i=1,nfo
	    Glob_NonlinParam(1:npt,nfru+i)=x_init((i-1)*npt+1:i*npt)
      enddo
      Glob_CurrEnergy=EnergyGA(nfrup1,cbs,.true.,ErrCode)
	  if (ErrCode/=0) then
        if (Glob_ProcID==0) write(*,*) 'Error in OptCycleG: energy cannot be computed'
		stop
	  endif
    endif

	if (CurrFunc>FuncEnd-NumOfFuncToShift) then
      LastIter=.true.
	else
      LastIter=.false.
	endif

    Glob_History(cbs)%Energy=Glob_CurrEnergy
	if (LastIter) then
      Glob_History(cbs)%InitFuncAtLastStep=0
	  Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
	else
      Glob_History(cbs)%InitFuncAtLastStep=CurrFunc
	endif

	if (Glob_ProcID==0) call SaveResults(Sort='yes')

  enddo !end cycle CurrCycle

  if (Glob_ProcID==0) then
    write(*,*)
	write(*,'(1x,a5,i4,a13)') 'Cycle',CurrCycle,' has finished'
  endif
  if (CurrCycle/=NumCycles) then
    Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    if (Glob_ProcID==0) write(*,'(1x,a47)',advance='no') 'Ordering basis functions and matrix elements...'
	call SortBasisFuncAndMatElem(fbn,cbs,FuncNumTemp,NonlinParamTemp,TempR,TempC)
    if (Glob_ProcID==0) write(*,*) 'done'
  endif

enddo !End of main optimization cycle 

if (Glob_ProcID==0) write(*,'(1x,a53)',advance='no') 'Final ordering basis functions and matrix elements...'
call SortBasisFuncAndMatElem(FuncBegin,cbs,FuncNumTemp,NonlinParamTemp,TempR,TempC)
call ReverseFuncOrder(FuncBegin,cbs)
call ReverseMatElemOrder(FuncBegin,cbs)
if (Glob_ProcID==0) write(*,*) 'done'

call StoreMatricesInSwapFile()

!deallocate workspace
deallocate(TempC)
deallocate(TempR)
deallocate(FuncNumTemp)
deallocate(NonlinParamTemp)
deallocate(grad) 
deallocate(x_best)
deallocate(x_init)
deallocate(x)

!deallocate arrays used by DRMNG
deallocate(V_init)
deallocate(V)
deallocate(D)

!deallocate workspace for SaveResults
deallocate(Glob_IntWorkArrForSaveResults)

!dellocate workspace for EnergyGB
deallocate(Glob_WkGR)

!Deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX)
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!Deallocate some global arrays
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_D)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) write (*,'(1x,a30)') 'Routine OptCycleG has finished'

end subroutine OptCycleG



subroutine OptCycleI(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                     NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold) 
!Subroutine OptCycleI performs an optimization of the basis set of K functions 
!by means of cyclic optimization of a set of NumOfFuncToOpt functions at a time.
!After a step is made, the subroutine takes another set of NumOfFuncToOpt 
!functions, whose numbers are shifted by NumFuncToShift (normally, one
!would probably want to use NumFuncToShift=NumOfFuncToOpt). This strategy
!is applied NumCycles times in a row to functions whose numbers range from 
!FuncBegin to FuncEnd. MaxEnergyEval is the limit of the energy evaluations 
!for each step. If it happens that at some step the basis functions being
!optimized become linerly dependent with other basis functions in the basis
!then such a change is rejected and this step is omitted. Pair linear 
!dependency is checked by comparing the absolute value of the overlap
!with OverlapThreshold. If OverlapThreshold<0 then no linear dependency
!check is performed. A similar check is performed for all linear parameters. 
!If adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions 
!are rejected and the procedure is repeted until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avaid this check it is enough to set LinCoeffThreshold to a value 
!equal or smaller than zero. Subroutines EnergyIA, EnergyIB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)     :: K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift
integer,intent(in)     :: NumCycles, MaxEnergyEval
real(dprec),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
integer      i,j,m,ip,AttemptToGetGoodOverlap,ii
integer      np,npt,nfo,nfa,nfru,nfrup1,nv,nvmax,nfco,fbn,cbs,nfs
integer      CurrCycle,CurrFunc,CurrFuncBegin
integer      q,nr,OptIterCounter
integer      ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
real(dprec)  Evalue,E_best
logical      IsSwapFileOK,ExitNeeded,LastIter
logical      IsOverlapBad,IsAnyLinCoeffBad
real(dprec)  t
real(dprec),allocatable,dimension(:,:) :: NonlinParamTemp
integer,allocatable,dimension(:)       :: FuncNumTemp
real(dprec),allocatable,dimension(:)   :: TempR
complex(dprec),allocatable,dimension(:):: TempC
real(dprec),allocatable,dimension(:)   :: x,grad,x_init,x_best
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V,V_init
integer,parameter    :: LIV=60
integer                 IV(LIV),IV_init(LIV)
integer                 LV
integer                 ALG

!Setting the values of some global variables
cbs=Glob_CurrBasisSize
Glob_GSEPSolutionMethod='I'
Glob_OverlapPenaltyAllowed=.false.
np=Glob_np
npt=Glob_npt
nvmax=NumOfFuncToOpt*npt
Glob_HSLeadDim=Glob_CurrBasisSize
Glob_HSBuffLen=Glob_CurrBasisSize*NumOfFuncToOpt
Glob_nfa=Glob_CurrBasisSize
nfco=FuncEnd-FuncBegin+1
fbn=FuncBegin+Glob_CurrBasisSize-FuncEnd

!Checking if cyclic optimization is already completed for this basis size
if ((Glob_History(cbs)%CyclesDone>=NumCycles).and. &
    (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd)) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine OptCycleI has started'
  write(*,'(1x,a14,i5)') 'Basis size is ',cbs
  write(*,'(1x,a38)') 'Cyclic optimization of basis functions'
  write(*,'(1x,i5,a9,i5,a21)') FuncBegin,'  through',FuncEnd,' is already completed'
  write(*,'(1x,a20)') 'Exiting OptCycleI...'
  write (*,'(1x,a30)') 'Routine OptCycleI has finished'
  return
endif

if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine OptCycleI has started'
  write(*,'(1x,a14,i5)') 'Basis size is ',cbs
  write(*,'(1x,a38)') 'Cyclic optimization of basis functions'
  write(*,'(1x,i5,a9,i5,a19)') FuncBegin,'  through',FuncEnd,'  will be performed'
  write(*,'(1x,a14,i7)') 'MaxEnergyEval=',MaxEnergyEval
endif

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_invD(cbs))
allocate(Glob_D(2*npt,NumOfFuncToOpt,cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GHEPIIS, which is called
!inside EnergyIA and EnergyIB
allocate(Glob_WorkForGHEPIIS(cbs))
allocate(Glob_LastEigvector(cbs))
Glob_LastEigvector(1:cbs)=ONE

!Allocate workspace for EnergyIB
allocate(Glob_WkGR(NumOfFuncToOpt*npt))

!Allocate workspace for SaveResults (it uses Sort='yes'
!option, which requires workspace)
allocate(Glob_IntWorkArrForSaveResults(cbs))

!Allocate arrays used by DRMNG
allocate(D(nvmax))
LV=71+nvmax*(nvmax+13)/2 + 1
allocate(V(LV))
allocate(V_init(LV))

!Allocate workspace
allocate(x(nvmax))
allocate(x_init(nvmax))
allocate(x_best(nvmax))
allocate(grad(nvmax))
allocate(NonlinParamTemp(npt,cbs-FuncBegin+1))
allocate(FuncNumTemp(cbs-FuncBegin+1))
allocate(TempR(cbs-FuncBegin+1))
allocate(TempC(cbs-FuncBegin+1))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again
	
!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
call DIVSET(ALG,IV_init,LIV,LV,V_init)
IV_init(17)=1000000 
IV_init(18)=1000000
IV_init(19)=0 !set summary print format
IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
V_init(31)=0.D0*ONE
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
IV_init(1)=12 !DIVSET has been called and some default values were changed

!Read swap file (if necessary) and distribute the data
call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Changing the order of the functions to be optimized and, if necessary, 
!the corresponding matrix elements to reverse
call ReverseFuncOrder(FuncBegin,FuncEnd)
if (IsSwapFileOK) call ReverseMatElemOrder(FuncBegin,FuncEnd)

!Shifting the set of basis functions to be optimized to
!the very end and, if necessary, doing the permutation of matrix elements
!to reflect this change.
call PermuteFunctions(FuncBegin,FuncEnd,FuncNumTemp,NonlinParamTemp)
if (IsSwapFileOK) call PermuteMatrixElements(FuncBegin,FuncEnd,TempR,TempC)

!If the last optimized function is greater than FuncEnd-1 or smaller than FuncBegin
!we need to change it so that the optimization begins from FuncBegin
if (Glob_History(cbs)%InitFuncAtLastStep<FuncBegin) Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
if (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd) then
  Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
  Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
endif

!We need to make an initial permutation to shift the functions that were 
!already optimized (if any) in the current optimization cycle
ip=Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
if (ip>0) then
  call PermuteFunctions(fbn,cbs-ip,FuncNumTemp,NonlinParamTemp)
  if (IsSwapFileOK) call PermuteMatrixElements(fbn,cbs-ip,TempR,TempC)   
endif

!Calculating the initial energy
call linalg_setparam(cbs)
if (IsSwapFileOK) then
  !Getting initial energy
  if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyIA(1,cbs,.false.,ErrCode)
else
  !Getting initial energy
  if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyIA(1,cbs,.true.,ErrCode)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in OptCycleI: initial energy cannot be computed'
  stop
endif

if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy


!Here comes main optimization cycle
do CurrCycle=Glob_History(cbs)%CyclesDone+1,NumCycles
  !Doing cycle number CurrCycle
  if (Glob_ProcID==0) then
    write(*,*)  
    write(*,'(1x,a5,i4,a10)') 'Cycle',CurrCycle,' has begun'
  endif
  CurrFuncBegin=Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
  q=FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
  OptIterCounter=0
  do CurrFunc=CurrFuncBegin,FuncEnd,NumOfFuncToShift
    !Note that CurrFunc counts functions as they were not shifted back
	!by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
	!initial input file.
    nfo=min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
	Glob_nfo=nfo	
    nv=nfo*npt
	nfru=cbs-nfo
	Glob_nfru=nfru
	nfrup1=nfru+1
	OptIterCounter=OptIterCounter+1
	Glob_InvItTempCounter1=0
    Glob_InvItTempCounter2=0	
    if (Glob_ProcID==0) then
      write(*,*)
	  if (nfo>1) then
	    write(*,'(1x,a21,i5,a10,i5)') 'Optimizing functions ',CurrFunc,'  through ',CurrFunc+nfo-1
      else
	    write(*,'(1x,a21,i5)') 'Optimizing function ',CurrFunc
	  endif
    endif   
	
    !Permute basis functions so that the last nfo functions, which have just
	!been optimized are replaced with those that will be optimized at next step	 
	nfs=min(NumOfFuncToShift,nfo)
    if (OptIterCounter/=1) then
      i=cbs-nfo-nfs+1
	  j=cbs-nfs
      call PermuteFunctions(i,j,FuncNumTemp,NonlinParamTemp)
	  call PermuteMatrixElements(i,j,TempR,TempC)
	  nr=NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)
	  if (q-nfo>=nr*2) then
        !normal permutation, there is sufficient number of functions left
		m=cbs-nfo
		j=m-nr
		i=j-nr+1       
      else
        !abnormal permutation, not sufficient number of functions left
		m=cbs-nfo
		j=m-nr
		i=cbs-q+1
	  endif
      call PermuteFunctions2(i,j,m,FuncNumTemp,NonlinParamTemp)
	  call PermuteMatrixElements2(i,j,m,TempR,TempC) 
	  !call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
	  Glob_CurrEnergy=EnergyIA(i,cbs,.false.,ErrCode)
      if (ErrCode/=0) then
        if (Glob_ProcID==0) write(*,*) 'Error in OptCycleI: energy cannot be computed after permuting basis functions'
        stop
      endif	  
    endif

    if (Glob_ProcID==0) then
	  if (Glob_AreParamPrintedInCycleOptX) then
        write (*,*) 'Nonlinear parameters before optimization:'
        do i=1,nfo
          write(*,'(1x,i5,a2,<npt>d23.16)') Glob_FuncNum(nfru+i),': ',Glob_NonlinParam(1:npt,nfru+i)
        enddo	
	  endif		    
    endif

    !Setting IV and V values as was in their initial copies
	IV(1:LIV)=IV_init(1:LIV)
	V(1:LV)=V_init(1:LV)
	do i=1,nfo
	  x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	  x_init((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
	enddo
	
	t=max(ONE/(cbs*cbs*sqrt(ONE*cbs)),10000*epsilon(Glob_CurrEnergy))
    do i=1,nfo
      !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
      do j=1,np
        !Make sure none of the D(i) will be zero or smaller than the threshold
        !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        D(npt*(i-1)+j)=t
      enddo
      !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
      D(npt*(i-1)+np+1:npt*i)=t
    enddo

	ExitNeeded=.false.
	NumOfFailures=0
    NumOfEnergyEval=0
    NumOfGradEval=0
    if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
    E_best=Glob_CurrEnergy
    x_best(1:nfo*npt)=x(1:nfo*npt)    

    do while (.not.(ExitNeeded))  
	  if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
      call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
	  select case (IV(1))
      case (1) !Only energy is needed
        call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	    do i=1,nfo
	      Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	    enddo
	    Evalue=EnergyIA(nfrup1,cbs,.true.,ErrCode)
        NumOfEnergyEval=NumOfEnergyEval+1
		if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
		  IV(2)=0
		else
          Glob_CurrEnergy=Evalue
          if (Evalue<E_best) then 
            E_best=Evalue
            x_best(1:nfo*npt)=x(1:nfo*npt)
          endif           
		endif
	  case (2) !Only gradient is needed
        call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	    do i=1,nfo
	      Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
	    enddo
        call EnergyIB(Evalue,grad,.true.,ErrCode)
        NumOfGradEval=NumOfGradEval+1
		if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
		  IV(2)=0
		else
          if (Evalue<E_best) then 
            E_best=Evalue
            x_best(1:nfo*npt)=x(1:nfo*npt)
          endif 			  
		endif	     
	  case (3:8) !Some kind of convergence has been reached
        ExitNeeded=.true.
      case (9:10) !Function evaluation limit has been reached.
	              !This is never supposed to happen because we
				  !count the number of function evaluations ourselves.
        ExitNeeded=.true.
	  endselect  
	  if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	    if (Glob_ProcID==0) then
          write(*,*) 'Error in OptCycleI: number of failures in energy or gradient'
		  write(*,*) 'calculations during the optimization of nonlinear parameters'
		  write(*,*) 'exceeded limit' 
		endif
	    stop
	  endif
	  if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
    enddo !while

    !Compute the energy and the linear coefficients at the best point found
    do i=1,nfo
	  Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
    enddo
    Glob_CurrEnergy=EnergyIAM(nfrup1,cbs,.true.,ErrCode)   
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
	    write(*,*) 'Warning in OptCycleI: failed to evaluate energy after optimization'
	    write(*,*) 'of nonlinear parameters. The values of the nonlinear parameters'
		write(*,*) 'will be left unchanged.'
	  endif 
    endif
 
 	!Checking if overlap is OK (only in case OverlapThreshold>ZERO)
    IsOverlapBad=.false.
    if (OverlapThreshold>ZERO) then
      ii=0
      do i=nfrup1,cbs
        do j=1,i-1
          if (abs(Glob_S(i,j))>OverlapThreshold) then
            ii=ii+1
			IsOverlapBad=.true.
	        if (Glob_ProcID==0) then
	          if (ii==1) then
                write(*,*) 'Warning: overlap of the following functions exceeds threshold'
		        write(*,*) 'Nonlinear parameters will be left unchanged'
		      endif
	          write(*,'(i5,a1,i5,i5,a7,d23.16,a1,d23.16,a1)') &
		           ii,':',i,j,'    S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),')'
            endif
	      endif
        enddo
      enddo
    endif
	!Checking if liniar coefficients are OK (only in case LinCoeffThreshold>ZERO)
	IsAnyLinCoeffBad=.false.
	if (LinCoeffThreshold>ZERO) then
      ii=0
	  do i=1,cbs
        if(abs(Glob_c(i))>LinCoeffThreshold) then
          ii=ii+1
		  IsAnyLinCoeffBad=.true.
		  if (Glob_ProcID==0) then
            if (ii==1) then
              write(*,*) 'Warning: absolute value of linear parameters of the'
			  write(*,*) 'following functions exceeds threshold'
		      write(*,*) 'Nonlinear parameters will be left unchanged'
			endif
			write(*,'(i5,a1,i5,a7,d23.16,a1,d23.16,a1)') &
			      ii,':',i,'    c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')'
		  endif
		endif
      enddo
	endif

    if ((Glob_ProcID==0).and.(ErrCode==0).and.(.not.IsOverlapBad).and.(.not.IsAnyLinCoeffBad)) then
      write (*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
      write (*,'(1x,a2,d23.16)') 'E=',Glob_CurrEnergy
      if (Glob_AreParamPrintedInCycleOptX) then
        write (*,*) 'Nonlinear parameters after optimization:'
        do i=1,nfo
          write(*,'(1x,i5,a2,<npt>d23.16)') Glob_FuncNum(nfru+i),': ',Glob_NonlinParam(1:npt,nfru+i)
        enddo
        write (*,*) 'Average number of iterations in GHEPIIS: ',(Glob_InvItTempCounter2*1.0)/Glob_InvItTempCounter1	  	    		    
      endif
    endif

    if ((ErrCode/=0).or.IsOverlapBad.or.IsAnyLinCoeffBad) then
	  !restore initial values of nonlinear parameters
      do i=1,nfo
	    Glob_NonlinParam(1:npt,nfru+i)=x_init((i-1)*npt+1:i*npt)
      enddo
      Glob_CurrEnergy=EnergyIA(nfrup1,cbs,.true.,ErrCode)
	  if (ErrCode/=0) then
        if (Glob_ProcID==0) write(*,*) 'Error in OptCycleI: energy cannot be computed'
		stop
	  endif
    endif

	if (CurrFunc>FuncEnd-NumOfFuncToShift) then
      LastIter=.true.
	else
      LastIter=.false.
	endif

    Glob_History(cbs)%Energy=Glob_CurrEnergy
	if (LastIter) then
      Glob_History(cbs)%InitFuncAtLastStep=0
	  Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
	else
      Glob_History(cbs)%InitFuncAtLastStep=CurrFunc
	endif

	if (Glob_ProcID==0) call SaveResults(Sort='yes')

  enddo !end cycle CurrCycle

  if (Glob_ProcID==0) then
    write(*,*)
	write(*,'(1x,a5,i4,a13)') 'Cycle',CurrCycle,' has finished'
  endif
  if (CurrCycle/=NumCycles) then
    Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    if (Glob_ProcID==0) write(*,'(1x,a47)',advance='no') 'Ordering basis functions and matrix elements...'
	call SortBasisFuncAndMatElem(fbn,cbs,FuncNumTemp,NonlinParamTemp,TempR,TempC)
	!call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
	Glob_CurrEnergy=EnergyIA(fbn,cbs,.false.,ErrCode)
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error in OptCycleI: energy cannot be computed after sorting basis functions'
      stop
    endif
    if (Glob_ProcID==0) write(*,*) 'done'
  endif

enddo !End of main optimization cycle 

if (Glob_ProcID==0) write(*,'(1x,a53)',advance='no') 'Final ordering basis functions and matrix elements...'
call SortBasisFuncAndMatElem(FuncBegin,cbs,FuncNumTemp,NonlinParamTemp,TempR,TempC)
call ReverseFuncOrder(FuncBegin,cbs)
call ReverseMatElemOrder(FuncBegin,cbs)
if (Glob_ProcID==0) write(*,*) 'done'

call StoreMatricesInSwapFile()

!deallocate workspace
deallocate(TempC)
deallocate(TempR)
deallocate(FuncNumTemp)
deallocate(NonlinParamTemp)
deallocate(grad) 
deallocate(x_best)
deallocate(x_init)
deallocate(x)

!deallocate arrays used by DRMNG
deallocate(V_init)
deallocate(V)
deallocate(D)

!deallocate workspace for SaveResults
deallocate(Glob_IntWorkArrForSaveResults)

!dellocate workspace for EnergyIB
deallocate(Glob_WkGR)

!Deallocate workspace for subroutine GHEPIIS
deallocate(Glob_LastEigvector)
deallocate(Glob_WorkForGHEPIIS)

!Deallocate some global arrays
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_D)
deallocate(Glob_invD)
deallocate(Glob_diagS)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) write (*,'(1x,a30)') 'Routine OptCycleI has finished'

end subroutine OptCycleI




subroutine FullOpt1G(InitFunc,FinalFunc,MaxEnergyEval,OverlapThreshold,MaxOverlapPenalty,  &
    DataSaveMinTimeInterv,HessianSaveMinTimeInterv,HessFileName)
!Subroutine FullOpt1G performs simultaneous optimization of 
!nonlinear parameters of basis functions whose number
!ranges from InitFunct to FinalFunc. Optimization routine
!used in FullOpt1G is DRMNG. Parameter MaxEnergyEval defines 
!the limit of the energy evaluations (just the energy, not
!the gradient) in the optimization process. The optimization uses penalty
!functions for the overlaps. If an overlap, Sij, exceeds the threshold, 
!OverlapThreshold, then a penaly is applied for this (the penalty function is 
!quadrativ and smooth so that the actual objective function is differentiable, 
!although the second derivative is not smooth). If the value of OverlapThreshold 
!is set equal or higher than 1.0 then no penalty will be applied. Parameter 
!MaxOverlapPenalty defines the magnitude of penalty. When MaxOverlapPenalty=1.0 
!and Sij=1.0 then the penalty due to this pair overlap will be 1.0. The routine
!saves the hessian in the file whose name is defined by 
!variable HessFileName (may include the full path if necessary).
!This save occurs periodically, with the interval approximately equal
!to the value of parameter HessianSaveMinTimeInterv. 
!Similarly, the input/output file (which contains basis set) is 
!updated periodically, with the interval approximately equal
!to the value of parameter DataSaveMinTimeInterv. 
!If there is not need to save the hessian at all, one can set 
!HessFileName=' ' or HessFileName='none' or HessFileName='NONE'.
!Any other name will be treated as a file name. It is important to 
!mention that the hessian itself and the corresponding file
!may require a lot of memory/storage (approximately k*k/2 real values, where
!k is the dimensionality of the problem; k=CurrBasisSize*Glob_npt).
!Onnce the routine begins optimization, it actually tries
!to locate file HessFileName and read the hessian from there. If this
!operation is succesful, this hessian is used to initialize
!routine DRMNG. Subroutines EnergyGA, EnergyGB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)      :: InitFunc,FinalFunc,MaxEnergyEval
real(dprec),intent(in)  :: OverlapThreshold,MaxOverlapPenalty,DataSaveMinTimeInterv,HessianSaveMinTimeInterv
character(Glob_FileNameLength),intent(in) :: HessFileName
!Local variables:
integer      i,j
integer      np,npt,nfo,nfa,nfru,nv
integer      OpenFileErr,ErrCode,NumOfEnergyEval,NumOfGradEval,NumOfFailures
integer      NumOfEnergyEvalDuringFullOpt_Init
integer      BlockSizeForZHEGVX
logical      IsSwapFileOK,ExitNeeded
logical      SaveHessian,IsHessFileOK,IsHessSaveSuccess
real(dprec)  Evalue,CurrentEnergy
real(dprec)  TimeOfLastSave, TimeOfLastHessSave
real(dprec)  t
integer      tas,InitFuncNew
complex(dprec)  MaxAbsOverlap,MinAbsOverlap
real(dprec)     AverageAbsOverlap
real(dprec),allocatable,dimension(:,:)   :: NonlinParamTemp
integer,allocatable,dimension(:)         :: FuncNumTemp
real(dprec),allocatable,dimension(:)     :: TempR
complex(dprec),allocatable,dimension(:)  :: TempC
real(dprec),allocatable,dimension(:)     :: x,grad 
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V
integer,parameter    :: LIV=60
integer                 IV(LIV)
integer                 LV
integer                 ALG
integer                 IVLMAT
!!==================================================== 
!!These variables are used when a finite difference gradient is computed               
!real(dprec)                                    deltax,Evalue1
!integer                                        m
!complex(dprec)                                 Skl_plus, Skl_minus
!real(dprec)                                    gar1,gai1,gar2,gai2,gnr,gni 
!!====================================================

if (OverlapThreshold>=ONE) then
  Glob_OverlapPenaltyAllowed=.false. 
else
  Glob_OverlapPenaltyAllowed=.true. 
  Glob_OverlapPenaltyThreshold2=OverlapThreshold*OverlapThreshold
  Glob_MaxOverlapPenalty=MaxOverlapPenalty
endif  
if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine FullOpt1G has started'
  write(*,'(1x,a68)') 'Simultaneous optimization of nonlinear parameters of basis functions'
  write(*,'(1x,i5,a9,i5,a19)') InitFunc,'  through',FinalFunc,'  will be attempted'
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a21,d23.16)') 'Overlap threshold is ',abs(OverlapThreshold)
    write(*,'(1x,a39,d23.16)') 'Max value of a pair overlap penalty is ',Glob_MaxOverlapPenalty
    write(*,'(1x,a68)') 'Warning! The energy value that will be shown during the optimization' 
    write(*,'(1x,a33)') 'may differ from the actual energy'    
  else
    write(*,'(1x,a41)') 'No constrains on overlaps will be imposed'
  endif  
endif

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_nfa=Glob_CurrBasisSize
Glob_nfo=FinalFunc-InitFunc+1
Glob_nfru=Glob_CurrBasisSize-Glob_nfo
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
nfa=Glob_nfa
nfo=Glob_nfo
nfru=Glob_nfru
nv=nfo*npt
InitFuncNew=nfa+InitFunc-FinalFunc
Glob_HSBuffLen=max(min(nfa*(nfa+1)/2,1000),30*nfa)

!Allocate some global arrays
allocate(Glob_H(nfa,nfa))
allocate(Glob_S(nfa,nfa))
allocate(Glob_diagH(nfa))
allocate(Glob_diagS(nfa))
allocate(Glob_D(2*npt,nfo,nfa))
allocate(Glob_c(nfa))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',nfa,nfa,nfa,nfa)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*nfa
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*nfa))
allocate(Glob_IWorkForZHEGVX(5*nfa)) 

!Allocate workspace for EnergyGB
allocate(Glob_WkGR(nfo*npt))

!Allocate arrays used by DRMNG
LV=71 + nv*(nv+13)/2 + 1
if (Glob_ProcID==0) then
  allocate(D(nv))
  allocate(V(LV))
endif

!Allocate workspace
allocate(x(nv))
allocate(grad(nv)) 

!Setting up a logical variable that determines
!whether the hessian should be saved from time to time
if ((HessFileName==' ').or.(HessFileName=='none').or. &
    (HessFileName=='NONE').or.(HessFileName=='None')) then
  SaveHessian=.false.
else
  SaveHessian=.true.
endif

call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Shifting basis functions InitFunc through FinalFunc to
!the very end and, if necessary, doing the permutation of 
!matrix elements to reflect this change in function order.
if (FinalFunc/=nfa) then
  !allocate space
  tas=FinalFunc-InitFunc+1
  allocate(NonlinParamTemp(1:npt,1:tas))
  allocate(FuncNumTemp(1:tas))
  allocate(TempR(1:tas))
  allocate(TempC(1:tas))
  call PermuteFunctions(InitFunc,FinalFunc,FuncNumTemp,NonlinParamTemp)
  if (IsSwapFileOK) call PermuteMatrixElements(InitFunc,FinalFunc,TempR,TempC)
  deallocate(TempC)
  deallocate(TempR)
  deallocate(FuncNumTemp)
  deallocate(NonlinParamTemp)
  !Allocate workspace for SaveResults (it uses Sort='yes'
  !option, which requires workspace)
  allocate(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
endif

!Calculating the initial energy
if (IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.false.,ErrCode)
else
  if (Glob_ProcID==0) write(*,'(1x,a59)',advance='no') 'Computing matrix elements and solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.true.,ErrCode)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in FullOpt1G: initial energy cannot be computed'
  stop
endif
if (Glob_ProcID==0) write(*,*) ' done'
if (Glob_ProcID==0) then
  call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a43,d23.16)') 'Initial energy (without overlap penalty)   ',Glob_CurrEnergy-Glob_TotalOverlapPenalty
    write(*,'(1x,a43,d23.16)') 'Overlap penalty                            ',Glob_TotalOverlapPenalty
    write(*,'(1x,a43,d23.16)') 'Initial energy (including overlap penalty) ',Glob_CurrEnergy      
  else
    write(*,'(1x,a43,d23.16)') 'Initial energy                             ',Glob_CurrEnergy
  endif
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
  write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap 
endif

Glob_TimeSinceStart=Glob_TimeSinceStart+DTIME(Glob_TmArray)
TimeOfLastSave=Glob_TimeSinceStart
TimeOfLastHessSave=Glob_TimeSinceStart

!Setting parameters for DRMNG
	
!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
if (Glob_ProcID==0) then
  call DIVSET(ALG,IV,LIV,LV,V)
  IV(17)=1000000; IV(18)=1000000
  IV(19)=-1 !set summary print format
  IV(20)=0; IV(22)=0; IV(23)=-1; IV(24)=0
  V(31)=0.D0*ONE
  !V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
  !VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
  !MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
  V(35)=Glob_MaxScStepAllowedInOpt
  IV(1)=12 !DIVSET has been called and some default values were changed
  nv=nfo*npt
endif
do i=1,nfo
    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,InitFuncNew+i-1)
enddo

if (Glob_ProcID==0) then
  !If SaveHessian=.true. then try to read the Hessian
  !from the file
  if (SaveHessian) then
    IVLMAT=IV(42)
    call ReadHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessFileOK)
    if (IsHessFileOK) IV(25)=0
  endif
  if ((.not.Glob_FullOptSaveD).or.(.not.IsHessFileOK).or.(.not.SaveHessian)) then
    !Set the scaling vector
    t=max(ONE/(nfa*nfa*sqrt(ONE*nfa)),10000*epsilon(Glob_CurrEnergy))
    do i=1,nfo
      !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
      do j=1,np
        !Make sure none of the D(i) will be zero or smaller than the threshold
        !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        D(npt*(i-1)+j)=t
      enddo
      !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
      D(npt*(i-1)+np+1:npt*i)=t
    enddo
  endif
endif

ExitNeeded=.false.
NumOfFailures=0
NumOfEnergyEval=0
NumOfGradEval=0
NumOfEnergyEvalDuringFullOpt_Init=Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt

do while (.not.(ExitNeeded))  
  if (Glob_ProcID==0) call DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
  call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  select case (IV(1))
  case (1) !Only energy is needed
    call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	do i=1,nfo
	  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
	enddo
	Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
    NumOfEnergyEval=NumOfEnergyEval+1
    if (ErrCode/=0) then
      NumOfFailures=NumOfFailures+1
	  IV(2)=1
	else
      CurrentEnergy=Evalue
	endif
    if (Glob_ProcID==0) then 
	  if (Evalue<Glob_CurrEnergy) then
        !Save data and Hessian if necessary
        Glob_TimeSinceStart=Glob_TimeSinceStart+DTIME(Glob_TmArray)
        if (Glob_TimeSinceStart-TimeOfLastSave>DataSaveMinTimeInterv) then
		  !Save the results if more than DataSaveMinTimeInterv seconds
		  !have passed since the last save
		  if (Glob_OverlapPenaltyAllowed) then
            Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
          else
            Glob_CurrEnergy=Evalue
          endif 
          !Changing history
          Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
          Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
             NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
          if (FinalFunc==nfa) then
            call SaveResults(Sort='no')
          else
            call SaveResults(Sort='yes')
		  endif
		  write(*,*) 'Data file has been updated'
          call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
          write(*,*) 'Some current statistics:'	
          if (Glob_OverlapPenaltyAllowed) then
            write(*,'(1x,a35,d23.16)') 'Energy (without overlap penalty)   ',Evalue-Glob_TotalOverlapPenalty
            write(*,'(1x,a35,d23.16)') 'Overlap penalty                    ',Glob_TotalOverlapPenalty
            write(*,'(1x,a35,d23.16)') 'Energy (including overlap penalty) ',Evalue        
          else
            write(*,'(1x,a35,d23.16)') 'Energy                             ',Evalue
          endif
          write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
          write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
          write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap 			  
          TimeOfLastSave=Glob_TimeSinceStart
        endif
        if ((Glob_TimeSinceStart-TimeOfLastHessSave>HessianSaveMinTimeInterv) &
		      .and.(SaveHessian)) then
		  !Save the Hessian if more than HessianSaveMinTimeInterv seconds
		  !have passed since the last save
          call SaveHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessSaveSuccess)
          TimeOfLastHessSave=Glob_TimeSinceStart
        endif
	  endif
    endif
  case (2) !Only gradient is needed
    call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	do i=1,nfo
	  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
	enddo
    call EnergyGB(Evalue,grad,.true.,ErrCode)
    NumOfGradEval=NumOfGradEval+1
    if (ErrCode/=0) then
      NumOfFailures=NumOfFailures+1
	  IV(2)=0
	endif
!	!===================================
!	!The lines below need be uncommented when finite
!    !difference gradient is shown 
!      if (Glob_ProcID==0) write(*,*) '     analytic gradient       finite diff gradient:'
!      do i=1,nfo
!        do j=1,npt
!          deltax=Glob_NonlinParam(j,i+nfru)*1.0D-6
!          Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)+deltax
!          Evalue1=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!          Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)-2*deltax
!          Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!          if (Glob_ProcID==0) write(*,'(1x,i3,1x,d23.16,1x,d23.16)') (i-1)*npt+j,grad((i-1)*npt+j),(Evalue1-Evalue)/(2*deltax)
!          Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)+deltax
!        enddo
!      enddo
!      if (Glob_ProcID==0) write(*,*)
!      if (Glob_ProcID==0) open(1,file='grad_S.txt',status='replace')
!      !output: analytic then numeric grad of Sij
!      i=6
!      j=1
!      !call SaveCMatrixSE(7,glob_S,'glob_S.txt')
!      if (Glob_ProcID==0) write(1,'(a5,d23.16,a1,d23.16,a1)') 'Sij=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),')'
!      if (Glob_ProcID==0) write(1,*)
!      if (Glob_ProcID==0) write(1,*) 'dSij/dParami'
!      do m=1,npt
!        deltax=Glob_NonlinParam(m,i)*1.0D-6
!        Glob_NonlinParam(m,i)=Glob_NonlinParam(m,i)+deltax
!        Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!        Skl_plus=Glob_S(i,j)
!        Glob_NonlinParam(m,i)=Glob_NonlinParam(m,i)-2*deltax
!        Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode) 
!        Skl_minus=Glob_S(i,j)
!        gar1=real(Glob_D(Glob_npt+m,i-Glob_nfru,j))
!        gai1=imag(Glob_D(Glob_npt+m,i-Glob_nfru,j))
!        gar2=real(Glob_D(Glob_npt+m,i-Glob_nfru,j)-ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(i,j)*sqrt(Glob_diagS(i))/sqrt(Glob_diagS(j)))
!        gai2=imag(Glob_D(Glob_npt+m,i-Glob_nfru,j)-ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(i,j)*sqrt(Glob_diagS(i))/sqrt(Glob_diagS(j)))
!        gnr=real(Skl_plus-Skl_minus)/(2.0D0*deltax)
!        gni=imag(Skl_plus-Skl_minus)/(2.0D0*deltax)        
!        if (Glob_ProcID==0) write(1,'(a1,d23.16,a1,d23.16,a4,d23.16,a1,d23.16,a4,d23.16,a1,d23.16,a1)') '(',gar1,',',gai1,')  (',gar2,',',gai2,')  (',gnr,',',gni,')'       
!        Glob_NonlinParam(m,i)=Glob_NonlinParam(m,i)+deltax         
!      enddo
!      write(1,*) 'dSij/dParamj'
!      do m=1,npt
!        deltax=Glob_NonlinParam(m,j)*1.0D-6
!        Glob_NonlinParam(m,j)=Glob_NonlinParam(m,j)+deltax
!        Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!        Skl_plus=Glob_S(i,j)
!        Glob_NonlinParam(m,j)=Glob_NonlinParam(m,j)-2*deltax
!        Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode) 
!        Skl_minus=Glob_S(i,j)
!        gar1=real(Glob_D(Glob_npt+m,j-Glob_nfru,i))
!        gai1=imag(Glob_D(Glob_npt+m,j-Glob_nfru,i))
!        gar2=real(Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*conjg(Glob_S(i,j))*sqrt(Glob_diagS(j))/sqrt(Glob_diagS(i)))
!        gai2=-imag(Glob_D(Glob_npt+m,j-Glob_nfru,i)-ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*conjg(Glob_S(i,j))*sqrt(Glob_diagS(j))/sqrt(Glob_diagS(i)))
!        gnr=real(Skl_plus-Skl_minus)/(2.0D0*deltax)
!        gni=imag(Skl_plus-Skl_minus)/(2.0D0*deltax)        
!        if (Glob_ProcID==0) write(1,'(a1,d23.16,a1,d23.16,a4,d23.16,a1,d23.16,a4,d23.16,a1,d23.16,a1)') '(',gar1,',',gai1,')  (',gar2,',',gai2,')  (',gnr,',',gni,')'               
!        Glob_NonlinParam(m,j)=Glob_NonlinParam(m,j)+deltax         
!      enddo      
!      if (Glob_ProcID==0) close(1)
!    !===================================	
  case (3:8) !Some kind of convergence has been reached
    ExitNeeded=.true.
  case (9:10) !Function evaluation limit has been reached.
	          !This never suppose to happen because we
			  !count the number of function evaluations ourselves.
    ExitNeeded=.true.
  endselect  
  if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	if (Glob_ProcID==0) then
      write(*,*) 'Error in FullOpt1G: number of failures in energy or gradient'
	  write(*,*) 'calculations during the optimization of nonlinear parameters'
	  write(*,*) 'exceeded limit' 
	endif
	stop
  endif
  if (NumOfEnergyEval>=MaxEnergyEval) then
    if (Glob_ProcID==0) then
      write(*,*) 'Warning in FullOpt1G: number of energy evaluations reached limit'
	  write(*,*) 'Optimization is terminated'
	endif
	ExitNeeded=.true.  
  endif
enddo

!Calculate the energy at the best point found
do i=1,nfo
  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
enddo
Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
if (ErrCode/=0) then
  if (Glob_ProcID==0) then
	write(*,*) 'Error in FullOpt1G: failed to evaluate energy after the optimization'
	write(*,*) 'of nonlinear parameters'
  endif      
  stop
endif
if (Glob_OverlapPenaltyAllowed) then
  Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
else
  Glob_CurrEnergy=Evalue
endif 

!Printing the number of energy/gradient evaluations
!and the energy after the optimization
if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
  write(*,*) 'Final energy and overlap statistics:'
  call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a35,d23.16)') 'Energy (without overlap penalty)   ',Evalue-Glob_TotalOverlapPenalty
    write(*,'(1x,a35,d23.16)') 'Overlap penalty                    ',Glob_TotalOverlapPenalty
    write(*,'(1x,a35,d23.16)') 'Energy (including overlap penalty) ',Evalue        
  else
    write(*,'(1x,a35,d23.16)') 'Energy                             ',Evalue
  endif
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
  write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap  
endif

!Adding data to history
Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
!Glob_History(Glob_CurrBasisSize)%CyclesDone=0
!Glob_History(Glob_CurrBasisSize)%InitFuncAtLastStep=0
Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
   NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval

!Shifting basis functions InitFunc through FinalFunc back from 
!the very end to the middle, where they initially were, and, if 
!necessary, doing the permutation of matrix elements to reflect 
!this change in function order.
if (FinalFunc/=nfa) then
  !allocate space
  tas=nfa-FinalFunc
  allocate(NonlinParamTemp(1:npt,1:tas))
  allocate(FuncNumTemp(1:tas))
  allocate(TempR(1:tas))
  allocate(TempC(1:tas))
  call PermuteFunctions(InitFunc,InitFunc+tas-1,FuncNumTemp,NonlinParamTemp)
  if (Glob_UseSwapFile) call PermuteMatrixElements(InitFunc,InitFunc+tas-1,TempR,TempC)
  deallocate(TempC)
  deallocate(TempR)
  deallocate(FuncNumTemp)
  deallocate(NonlinParamTemp)
  !deallocate workspace for SaveResults
  deallocate(Glob_IntWorkArrForSaveResults)
endif

!saving results
if (Glob_ProcID==0) call SaveResults(Sort='no')

call StoreMatricesInSwapFile()

!deallocate workspace
deallocate(grad) 
deallocate(x)

!deallocate arrays used by DRMNG
if (Glob_ProcID==0) then
  deallocate(V)
  deallocate(D)
endif
!dellocate workspace for EnergyGB
deallocate(Glob_WkGR)

!Deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX)
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!Deallocate some global arrays
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_D)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) write (*,'(1x,a30)') 'Routine FullOpt1G has finished'

end subroutine FullOpt1G



subroutine FullOpt1I(InitFunc,FinalFunc,MaxEnergyEval,OverlapThreshold,MaxOverlapPenalty,  &
    DataSaveMinTimeInterv,HessianSaveMinTimeInterv,HessFileName)
!Subroutine FullOpt1I performs simultaneous optimization of 
!nonlinear parameters of basis functions whose number
!ranges from InitFunct to FinalFunc. Optimization routine
!used in FullOpt1I is DRMNG. Parameter MaxEnergyEval defines 
!the limit of the energy evaluations (just the energy, not
!the gradient) in the optimization process. The optimization uses penalty
!functions for the overlaps. If an overlap, Sij, exceeds the threshold, 
!OverlapThreshold, then a penaly is applied for this (the penalty function is 
!quadrativ and smooth so that the actual objective function is differentiable, 
!although the second derivative is not smooth). If the value of OverlapThreshold 
!is set equal or higher than 1.0 then no penalty will be applied. Parameter 
!MaxOverlapPenalty defines the magnitude of penalty. When MaxOverlapPenalty=1.0 
!and Sij=1.0 then the penalty due to this pair overlap will be 1.0. The routine
!saves the hessian in the file whose name is defined by 
!variable HessFileName (may include the full path if necessary).
!This save occurs periodically, with the interval approximately equal
!to the value of parameter HessianSaveMinTimeInterv. 
!Similarly, the input/output file (which contains basis set) is 
!updated periodically, with the interval approximately equal
!to the value of parameter DataSaveMinTimeInterv. 
!If there is not need to save the hessian at all, one can set 
!HessFileName=' ' or HessFileName='none' or HessFileName='NONE'.
!Any other name will be treated as a file name. It is important to 
!mention that the hessian itself and the corresponding file
!may require a lot of memory/storage (approximately k*k/2 real values, where
!k is the dimensionality of the problem; k=CurrBasisSize*Glob_npt).
!Onnce the routine begins optimization, it actually tries
!to locate file HessFileName and read the hessian from there. If this
!operation is succesful, this hessian is used to initialize
!routine DRMNG. Subroutines EnergyIA, EnergyIB are
!called to evaluate the energy and its gradient.

!Arguments:
integer,intent(in)      :: InitFunc,FinalFunc,MaxEnergyEval
real(dprec),intent(in)  :: OverlapThreshold,MaxOverlapPenalty,DataSaveMinTimeInterv,HessianSaveMinTimeInterv
character(Glob_FileNameLength),intent(in) :: HessFileName
!Local variables:
integer      i,j
integer      np,npt,nfo,nfa,nfru,nv
integer      OpenFileErr,ErrCode,NumOfEnergyEval,NumOfGradEval,NumOfFailures
integer      NumOfEnergyEvalDuringFullOpt_Init
logical      IsSwapFileOK,ExitNeeded
logical      SaveHessian,IsHessFileOK,IsHessSaveSuccess
real(dprec)  Evalue,CurrentEnergy
real(dprec)  TimeOfLastSave, TimeOfLastHessSave
real(dprec)  t
integer      tas,InitFuncNew
complex(dprec)  MaxAbsOverlap,MinAbsOverlap
real(dprec)     AverageAbsOverlap
real(dprec),allocatable,dimension(:,:)   :: NonlinParamTemp
integer,allocatable,dimension(:)         :: FuncNumTemp
real(dprec),allocatable,dimension(:)     :: TempR
complex(dprec),allocatable,dimension(:)  :: TempC
real(dprec),allocatable,dimension(:)     :: x,grad 
!Arrays used by DRMNG
real(dprec),allocatable,dimension(:)     :: D,V
integer,parameter    :: LIV=60
integer                 IV(LIV)
integer                 LV
integer                 ALG
integer                 IVLMAT

if (OverlapThreshold>=ONE) then
  Glob_OverlapPenaltyAllowed=.false. 
else
  Glob_OverlapPenaltyAllowed=.true. 
  Glob_OverlapPenaltyThreshold2=OverlapThreshold*OverlapThreshold
  Glob_MaxOverlapPenalty=MaxOverlapPenalty
endif  
if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a29)') 'Routine FullOpt1I has started'
  write(*,'(1x,a68)') 'Simultaneous optimization of nonlinear parameters of basis functions'
  write(*,'(1x,i5,a9,i5,a19)') InitFunc,'  through',FinalFunc,'  will be attempted'
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a21)') 'Overlap threshold is ',abs(OverlapThreshold)
    write(*,'(1x,a39)') 'Max value of a pair overlap penalty is ',Glob_MaxOverlapPenalty
    write(*,'(1x,a68)') 'Warning! The energy value that will be shown during the optimization' 
    write(*,'(1x,a33)') 'may differ from the actual energy'    
  else
    write(*,'(1x,a41)') 'No constrains on overlaps will be imposed'
  endif  
endif

!Setting the values of some global variables
Glob_GSEPSolutionMethod='I'
Glob_nfa=Glob_CurrBasisSize
Glob_nfo=FinalFunc-InitFunc+1
Glob_nfru=Glob_CurrBasisSize-Glob_nfo
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
nfa=Glob_nfa
nfo=Glob_nfo
nfru=Glob_nfru
nv=nfo*npt
InitFuncNew=nfa+InitFunc-FinalFunc
Glob_HSBuffLen=max(min(nfa*(nfa+1)/2,1000),30*nfa)

!Allocate some global arrays
allocate(Glob_H(nfa,nfa))
allocate(Glob_S(nfa,nfa))
allocate(Glob_diagS(nfa))
allocate(Glob_invD(nfa))
allocate(Glob_D(2*npt,nfo,nfa))
allocate(Glob_c(nfa))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))
allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GHEPIIS, which is called
!inside EnergyIA and EnergyIB
allocate(Glob_WorkForGHEPIIS(nfa))
allocate(Glob_LastEigvector(nfa))
Glob_LastEigvector(1:nfa)=ONE

!Allocate workspace for EnergyIB
allocate(Glob_WkGR(nfo*npt))

!Allocate arrays used by DRMNG
LV=71 + nv*(nv+13)/2 + 1
if (Glob_ProcID==0) then
  allocate(D(nv))
  allocate(V(LV))
endif

!Allocate workspace
allocate(x(nv))
allocate(grad(nv)) 

!Setting up a logical variable that determines
!whether the hessian should be saved from time to time
if ((HessFileName==' ').or.(HessFileName=='none').or. &
    (HessFileName=='NONE').or.(HessFileName=='None')) then
  SaveHessian=.false.
else
  SaveHessian=.true.
endif

call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Shifting basis functions InitFunc through FinalFunc to
!the very end and, if necessary, doing the permutation of 
!matrix elements to reflect this change in function order.
if (FinalFunc/=nfa) then
  !allocate space
  tas=FinalFunc-InitFunc+1
  allocate(NonlinParamTemp(1:npt,1:tas))
  allocate(FuncNumTemp(1:tas))
  allocate(TempR(1:tas))
  allocate(TempC(1:tas))
  call PermuteFunctions(InitFunc,FinalFunc,FuncNumTemp,NonlinParamTemp)
  if (IsSwapFileOK) call PermuteMatrixElements(InitFunc,FinalFunc,TempR,TempC)
  deallocate(TempC)
  deallocate(TempR)
  deallocate(FuncNumTemp)
  deallocate(NonlinParamTemp)
  !Allocate workspace for SaveResults (it uses Sort='yes'
  !option, which requires workspace)
  allocate(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
endif

!Calculating the initial energy
call linalg_setparam(Glob_CurrBasisSize)
if (IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.false.,ErrCode)
else
  if (Glob_ProcID==0) write(*,'(1x,a59)',advance='no') 'Computing matrix elements and solving eigenvalue problem...'
  Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.true.,ErrCode)
endif
if (ErrCode/=0) then
  if (Glob_ProcID==0) write(*,*) 'Error in FullOpt1I: initial energy cannot be computed'
  stop
endif
if (Glob_ProcID==0) write(*,*) ' done'
if (Glob_ProcID==0) then
  call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a43,d23.16)') 'Initial energy (without overlap penalty)   ',Glob_CurrEnergy-Glob_TotalOverlapPenalty
    write(*,'(1x,a43,d23.16)') 'Overlap penalty                            ',Glob_TotalOverlapPenalty
    write(*,'(1x,a43,d23.16)') 'Initial energy (including overlap penalty) ',Glob_CurrEnergy      
  else
    write(*,'(1x,a43,d23.16)') 'Initial energy                             ',Glob_CurrEnergy
  endif
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
  write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap 
endif

Glob_TimeSinceStart=Glob_TimeSinceStart+DTIME(Glob_TmArray)
TimeOfLastSave=Glob_TimeSinceStart
TimeOfLastHessSave=Glob_TimeSinceStart

!Setting parameters for DRMNG
	
!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
ALG=2
if (Glob_ProcID==0) then
  call DIVSET(ALG,IV,LIV,LV,V)
  IV(17)=1000000; IV(18)=1000000
  IV(19)=-1 !set summary print format
  IV(20)=0; IV(22)=0; IV(23)=-1; IV(24)=0
  V(31)=0.D0*ONE
  !V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
  !VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
  !MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
  V(35)=Glob_MaxScStepAllowedInOpt
  IV(1)=12 !DIVSET has been called and some default values were changed
  nv=nfo*npt
endif
do i=1,nfo
    x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,InitFuncNew+i-1)
enddo

if (Glob_ProcID==0) then
  !If SaveHessian=.true. then try to read the Hessian
  !from the file
  if (SaveHessian) then
    IVLMAT=IV(42)
    call ReadHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessFileOK)
    if (IsHessFileOK) IV(25)=0
  endif
  if ((.not.Glob_FullOptSaveD).or.(.not.IsHessFileOK).or.(.not.SaveHessian)) then
    !Set the scaling vector
    t=max(ONE/(nfa*nfa*sqrt(ONE*nfa)),10000*epsilon(Glob_CurrEnergy))
    do i=1,nfo
      !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
      do j=1,np
        !Make sure none of the D(i) will be zero or smaller than the threshold
        !D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        D(npt*(i-1)+j)=t
      enddo
      !D(npt*(i-1)+np+1:npt*i)=1.0*ONE
      D(npt*(i-1)+np+1:npt*i)=t
    enddo
  endif
endif

ExitNeeded=.false.
NumOfFailures=0
NumOfEnergyEval=0
NumOfGradEval=0
NumOfEnergyEvalDuringFullOpt_Init=Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt

do while (.not.(ExitNeeded))  
  if (Glob_ProcID==0) call DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
  call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  select case (IV(1))
  case (1) !Only energy is needed
    call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	do i=1,nfo
	  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
	enddo
	Evalue=EnergyIA(InitFuncNew,nfa,.true.,ErrCode)
    NumOfEnergyEval=NumOfEnergyEval+1
    if (ErrCode/=0) then
      NumOfFailures=NumOfFailures+1
	  IV(2)=1
	else
      CurrentEnergy=Evalue
	endif
    if (Glob_ProcID==0) then 
	  if (Evalue<Glob_CurrEnergy) then
        !Save data and Hessian if necessary
        Glob_TimeSinceStart=Glob_TimeSinceStart+DTIME(Glob_TmArray)
        if (Glob_TimeSinceStart-TimeOfLastSave>DataSaveMinTimeInterv) then
		  !Save the results if more than DataSaveMinTimeInterv seconds
		  !have passed since the last save
		  if (Glob_OverlapPenaltyAllowed) then
            Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
          else
            Glob_CurrEnergy=Evalue
          endif 
          !Changing history
          Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
          Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
             NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
          if (FinalFunc==nfa) then
            call SaveResults(Sort='no')
          else
            call SaveResults(Sort='yes')
		  endif
		  write(*,*) 'Data file has been updated'
          call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
          write(*,*) 'Some current statistics:'			
          if (Glob_OverlapPenaltyAllowed) then
            write(*,'(1x,a35,d23.16)') 'Energy (without overlap penalty)   ',Evalue-Glob_TotalOverlapPenalty
            write(*,'(1x,a35,d23.16)') 'Overlap penalty                    ',Glob_TotalOverlapPenalty
            write(*,'(1x,a35,d23.16)') 'Energy (including overlap penalty) ',Evalue        
          else
            write(*,'(1x,a35,d23.16)') 'Energy                             ',Evalue
          endif
          write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
          write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
          write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap 	          
          TimeOfLastSave=Glob_TimeSinceStart
        endif
        if ((Glob_TimeSinceStart-TimeOfLastHessSave>HessianSaveMinTimeInterv) &
		      .and.(SaveHessian)) then
		  !Save the Hessian if more than HessianSaveMinTimeInterv seconds
		  !have passed since the last save
          call SaveHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessSaveSuccess)
          TimeOfLastHessSave=Glob_TimeSinceStart
        endif
	  endif
    endif
  case (2) !Only gradient is needed
    call MPI_BCAST(x,nv,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
 	do i=1,nfo
	  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
	enddo
    call EnergyIB(Evalue,grad,.true.,ErrCode)
    NumOfGradEval=NumOfGradEval+1
    if (ErrCode/=0) then
      NumOfFailures=NumOfFailures+1
	  IV(2)=0
	endif
  case (3:8) !Some kind of convergence has been reached
    ExitNeeded=.true.
  case (9:10) !Function evaluation limit has been reached.
	          !This never suppose to happen because we
			  !count the number of function evaluations ourselves.
    ExitNeeded=.true.
  endselect  
  if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
	if (Glob_ProcID==0) then
      write(*,*) 'Error in FullOpt1I: number of failures in energy or gradient'
	  write(*,*) 'calculations during the optimization of nonlinear parameters'
	  write(*,*) 'exceeded limit' 
	endif
	stop
  endif
  if (NumOfEnergyEval>=MaxEnergyEval) then
    if (Glob_ProcID==0) then
      write(*,*) 'Warning in FullOpt1I: number of energy evaluations reached limit'
	  write(*,*) 'Optimization is terminated'
	endif
	ExitNeeded=.true.  
  endif
enddo

!Calculate the energy at the best point found
do i=1,nfo
  Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
enddo
Evalue=EnergyIA(InitFuncNew,nfa,.true.,ErrCode)
if (ErrCode/=0) then
  if (Glob_ProcID==0) then
	write(*,*) 'Error in FullOpt1I: failed to evaluate energy after the optimization'
	write(*,*) 'of nonlinear parameters'
  endif      
  stop
endif
if (Glob_OverlapPenaltyAllowed) then
  Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
else
  Glob_CurrEnergy=Evalue
endif 

!Printing the number of energy/gradient evaluations
!and the energy after the optimization
if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) 'Number of energy/gradient evaluations',NumOfEnergyEval,NumOfGradEval
  write(*,*) 'Final energy and overlap statistics:'
  call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
  if (Glob_OverlapPenaltyAllowed) then
    write(*,'(1x,a35,d23.16)') 'Energy (without overlap penalty)   ',Evalue-Glob_TotalOverlapPenalty
    write(*,'(1x,a35,d23.16)') 'Overlap penalty                    ',Glob_TotalOverlapPenalty
    write(*,'(1x,a35,d23.16)') 'Energy (including overlap penalty) ',Evalue        
  else
    write(*,'(1x,a35,d23.16)') 'Energy                             ',Evalue
  endif
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Maximal overlap  (',real(MaxAbsOverlap),',',imag(MaxAbsOverlap),')  |S|=',abs(MaxAbsOverlap)
  write(*,'(1x,a18,d16.9,a1,d16.9,a7,f16.9)') &
                'Minimal overlap  (',real(MinAbsOverlap),',',imag(MinAbsOverlap),')  |S|=',abs(MinAbsOverlap)
  write(*,'(1x,a29,f16.9)') 'Average abs value of overlap ',AverageAbsOverlap    
endif

!Adding data to history
Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
!Glob_History(Glob_CurrBasisSize)%CyclesDone=0
!Glob_History(Glob_CurrBasisSize)%InitFuncAtLastStep=0
Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
   NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval

!Shifting basis functions InitFunc through FinalFunc back from 
!the very end to the middle, where they initially were, and, if 
!necessary, doing the permutation of matrix elements to reflect 
!this change in function order.
if (FinalFunc/=nfa) then
  !allocate space
  tas=nfa-FinalFunc
  allocate(NonlinParamTemp(1:npt,1:tas))
  allocate(FuncNumTemp(1:tas))
  allocate(TempR(1:tas))
  allocate(TempC(1:tas))
  call PermuteFunctions(InitFunc,InitFunc+tas-1,FuncNumTemp,NonlinParamTemp)
  if (Glob_UseSwapFile) call PermuteMatrixElements(InitFunc,InitFunc+tas-1,TempR,TempC)
  deallocate(TempC)
  deallocate(TempR)
  deallocate(FuncNumTemp)
  deallocate(NonlinParamTemp)
  !deallocate workspace for SaveResults
  deallocate(Glob_IntWorkArrForSaveResults)
endif

!saving results
if (Glob_ProcID==0) call SaveResults(Sort='no')

call StoreMatricesInSwapFile()

!deallocate workspace
deallocate(grad) 
deallocate(x)

!deallocate arrays used by DRMNG
if (Glob_ProcID==0) then
  deallocate(V)
  deallocate(D)
endif
!dellocate workspace for EnergyIB
deallocate(Glob_WkGR)

!Deallocate workspace for GHEPSII
deallocate(Glob_LastEigvector)
deallocate(Glob_WorkForGHEPIIS)

!Deallocate some global arrays
deallocate(Glob_DlBuff2)
deallocate(Glob_DlBuff1)
deallocate(Glob_DkBuff2)
deallocate(Glob_DkBuff1)
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_D)
deallocate(Glob_invD)
deallocate(Glob_diagS)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) write (*,'(1x,a30)') 'Routine FullOpt1I has finished'

end subroutine FullOpt1I




subroutine EliminateLittleContribFunc(LinCoeffThreshold,FileName,PrintInfoSpec)
!Subroutine EliminateLittleContribFunc eliminates basis
!functions whose contribution to the energy is small. More
!precisely, it eliminates functions that have linear coefficients
!whose absolute values are smaller than LinCoeffThreshold. It is
!important to note that this subroutine uses coeffcients in
!front of normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1).
!The result is stored in file whose name is defined by parameter
!FileName. After saving the results the subroutine terminates
!the program execution.
!Parameter PrintInfoSpec specifies what information should be
!shown about linear coefficients:
! PrintInfoSpec=0,1 : the subroutine does not print any info
!                   about linear coefficients.
! PrintInfoSpec=2   : the subroutine outputs linear coefficients
!                     of all functions

!Arguments:
real(dprec),intent(in)   :: LinCoeffThreshold
character(Glob_FileNameLength),intent(in) :: FileName
integer,intent(in)       :: PrintInfoSpec

!Local variables:
integer      i,j
integer      np,npt,cbs
integer      OpenFileErr,ErrorCode,IFAIL
logical      IsSwapFileOK
integer      BlockSizeForZHEGVX
integer      NumOfEigvalsFound
real(dprec)  Evalue, EVs(1)
complex(dprec) Min_c,Max_c
real(dprec)  Aver_c
real(dprec),allocatable,dimension(:,:)   :: NonlinParamTemp
character(Glob_FileNameLength)          :: ch_temp



if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a46)') 'Routine EliminateLittleContribFunc has started'
  write(*,'(1x,a38,d23.16)') 'Linear coefficient threshold value is ',LinCoeffThreshold
endif 

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
cbs=Glob_CurrBasisSize

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs)) 

!Reading data from swap file

call ReadSwapFileAndDistributeData(IsSwapFileOK)

if (IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
else
  if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...' 
  call ComputeMatElem(1,cbs)
  if (Glob_ProcID==0) write(*,*) ' done'
  if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in EliminateLittleContribFunc: initial energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*) 'Basis size before elimination ',cbs
  write(*,*) 'Energy before elimination     ',Evalue
endif

allocate(NonlinParamTemp(1:npt,cbs))

if ((PrintInfoSpec>1).and.(Glob_ProcID==0)) then
  write(*,*) 'List of all linear coefficients:'
  do i=1,cbs
	write(*,'(i5,2x,a3,d23.16,a1,d23.16,a6,f16.9)') i,'c=(',real(Glob_c(i)),',',imag(Glob_c(i)),') |c|=',abs(Glob_c(i))
  enddo
  write(*,*)
endif

Min_c=Glob_c(1)
Max_c=Glob_c(1)
Aver_c=ZERO
j=0
do i=1,cbs
  if (abs(Glob_c(i))>=LinCoeffThreshold) then
    NonlinParamTemp(1:npt,i-j)=Glob_NonlinParam(1:npt,i)
  else
    if ((j==0).and.(Glob_ProcID==0)) write(*,*) 'Little contributing function list:'
    j=j+1
	write(*,'(i4,a1,i6,2x,a3,d16.9,a1,d16.9,a7,f19.12)') j,':',i,'c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')  |c|=',abs(Glob_c(i))
  endif
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Aver_c=Aver_c+abs(Glob_c(i))
enddo
Aver_c=Aver_c/cbs

if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a24,d16.9,a1,d16.9,a7,d16.9)') &
    'Smallest linear coeff. (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a24,d16.9,a1,d16.9,a7,d16.9)') &
    'Largest linear coeff.  (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a40,d23.16)') 'Average absolute value of linear coeff. ',Aver_c
  write(*,*)
endif

if (j==0) then
  if (Glob_ProcID==0) then
    write(*,*) 'There are no functions with the contribution' 
	write(*,*) 'smaller than ',LinCoeffThreshold
    write(*,*) 'No file have been written. Program will now stop'
  endif
  stop
endif

if (Glob_ProcID==0) then
  write(*,*) 'Basis size before elimination ',cbs
  write(*,*) 'Energy before elimination     ',Evalue
endif

Glob_CurrBasisSize=cbs-j
cbs=Glob_CurrBasisSize
Glob_NonlinParam(1:npt,1:cbs)=NonlinParamTemp(1:npt,1:cbs)

if (Glob_ProcID==0) then
  write(*,*) 'Computing matrix elements and solving eigenvalue problem using the'
  write(*,*) 'basis where little contributing functions have been eliminated...'
endif
call ComputeMatElem(1,cbs)

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in EliminateLittleContribFunc: energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) 'Basis size after elimination  ',cbs
  write(*,*) 'Energy after elimination      ',Evalue
endif

do i=1,cbs
  Glob_History(i)%Energy=ZERO
  Glob_History(i)%CyclesDone=0
  Glob_History(i)%InitFuncAtLastStep=0
  Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
enddo
Glob_History(cbs)%Energy=Evalue

Glob_CurrEnergy=Evalue
Glob_LastEigvalTol=1.0D30
Glob_BestEigvalTol=1.0D30
Glob_WorstEigvalTol=1.0D-30

ch_temp=Glob_DataFileName
Glob_DataFileName=FileName
if (Glob_ProcID==0) call SaveResults(Sort='no')
Glob_DataFileName=ch_temp

deallocate(NonlinParamTemp)

!deallocate global arrays
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

!deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX) 
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

if (Glob_ProcID==0) then
  i=len_trim(FileName)
  write(*,*) 'New basis has been saved in file ',FileName(1:i)
  write(*,*) 'Program will now stop'
endif

stop

end subroutine EliminateLittleContribFunc



subroutine EliminateLinDepFunc1(LinDepThreshold,FileName,PrintInfoSpec)
!Subroutine EliminateLinDepFunc1 eliminates linearly dependent
!functions. It checks for pair linear dependency only. It 
!removes from the basis those functions whose overlap (absolute value)
!with any of other basis function of smaller number is greater
!than LinDepThreshold. It is important to note that this subroutine 
!uses normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1).
!The result is stored in file whose name is defined by parameter
!FileName. After saving the results the program is terminated.
!Parameter PrintInfoSpec specifies what information should be
!shown during linear dependency check:
! PrintInfoSpec=0 : the subroutine does not print any info
!                   about basis functions that are linearly dependent.
! PrintInfoSpec=1 : the subroutine prints the values of the overlap of 
!                   linerly dependent functions.
! PrintInfoSpec=2 : in addition to the latter case it also prints the
!                   nonlinear parameters of linearly dependent functions.

!Arguments:
real(dprec),intent(in)   :: LinDepThreshold
character(Glob_FileNameLength),intent(in) :: FileName
integer,intent(in)       :: PrintInfoSpec

!Local variables:
integer        i,j,k
integer        np,npt,cbs
integer        OpenFileErr,ErrorCode,IFAIL
logical        IsSwapFileOK
integer        BlockSizeForZHEGVX
integer        NumOfEigvalsFound
real(dprec)    Evalue, EVs(1)
complex(dprec) MaxOverlap,MinOverlap
real(dprec)    AverOverlap
complex(dprec) Min_c,Max_c
real(dprec)    Average_c
integer,allocatable,dimension(:)         :: MaskArray
character(Glob_FileNameLength)          :: ch_temp



if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a40)') 'Routine EliminateLinDepFunc1 has started'
  write(*,'(1x,a27,d23.16)') 'Overlap threshold value is ',LinDepThreshold
endif 

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
cbs=Glob_CurrBasisSize

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs)) 

!Allocate local workspace
allocate(MaskArray(1:cbs))

!Reading data from swap file

call ReadSwapFileAndDistributeData(IsSwapFileOK)

if (.not.IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
  call ComputeMatElem(1,cbs)
  if (Glob_ProcID==0) write(*,*) ' done'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'


if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in EliminateLinDepFunc1: initial energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*) 'Pair linear dependency check:'
  write(*,*)
endif
do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

!Check overlap
MaskArray(1:cbs)=0
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
    if (abs(Glob_S(j,i))>LinDepThreshold) then
       MaskArray(j)=MaskArray(j)+1
	   k=k+1
	   if (Glob_ProcID==0) then
	     if (PrintInfoSpec>0) write(*,'(i5,a1,i5,i5,a6,d16.9,a1,d16.9,a6,f11.7)') &
		       k,':',i,j,'   S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),') |S|=',abs(Glob_S(i,j))
         if (PrintInfoSpec==2) then
           write(*,'(i6,1x,<Glob_npt>(1x,d23.16))') i,Glob_NonlinParam(1:Glob_npt,i)
		   write(*,'(6x,a3,d16.9,a1,d16.9,a9,f16.9)') &
		     'c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')    |c|=',abs(Glob_c(i))
           write(*,'(i6,1x,<Glob_npt>(1x,d23.16))') j,Glob_NonlinParam(1:Glob_npt,j)
		   write(*,'(6x,a3,d16.9,a1,d16.9,a9,f16.9)')  &
		     'c=(',real(Glob_c(j)),',',imag(Glob_c(j)),')    |c|=',abs(Glob_c(j))
		 endif
         write(*,*)
       endif
	endif
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients:
Min_c=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) '========== Summary before elimination: =========='
endif

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
endif

if (Glob_ProcID==0) then
  write(*,*) 'Basis size before elimination ',cbs
  write(*,*) 'Energy before elimination     ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif

if (k==0) then
  if (Glob_ProcID==0) then
    write(*,*) 'No linearly dependent functions have been found' 
    write(*,*) 'No file have been written. Program will now stop'
  endif
  stop
endif

j=0
i=1
do while (i+j<=cbs)
  if (MaskArray(i+j)>0) then
    j=j+1
  else
    Glob_NonlinParam(1:npt,i)=Glob_NonlinParam(1:npt,i+j)
	i=i+1
  endif
enddo

Glob_CurrBasisSize=cbs-k
cbs=Glob_CurrBasisSize

if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
call ComputeMatElem(1,cbs)
if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in EliminateLinDepFunc1: energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*)
  write(*,*) '========== Summary after elimination: ==========='
endif

!Check overlap
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients:
Min_c=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
  write(*,*) 'Basis size after elimination  ',cbs
  write(*,*) 'Energy after elimination      ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif

do i=1,cbs
  Glob_History(i)%Energy=ZERO
  Glob_History(i)%CyclesDone=0
  Glob_History(i)%InitFuncAtLastStep=0
  Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
enddo
Glob_History(cbs)%Energy=Evalue

Glob_CurrEnergy=Evalue
Glob_LastEigvalTol=1.0D30
Glob_BestEigvalTol=1.0D30
Glob_WorstEigvalTol=1.0D-30

ch_temp=Glob_DataFileName
Glob_DataFileName=FileName
if (Glob_ProcID==0) call SaveResults(Sort='no')
Glob_DataFileName=ch_temp

!dellocate local workspace
deallocate(MaskArray)

!deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX) 
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!deallocate global arrays
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) then
  i=len_trim(FileName)
  write(*,*) 'New basis has been saved in file ',FileName(1:i)
  write(*,*) 'Program will now stop'
endif

stop

end subroutine EliminateLinDepFunc1


subroutine SeparateLinDepFunc1(LinDepThreshold,SeparationParam,FileName,PrintInfoSpec)
!Subroutine SeparateLinDepFunc1 does exactly the same thing as 
!subroutine EliminateLinDepFunc1 does, but without throwing away
!linearly dependent functions. Instead, it changes the parameters
!of such functions randomly (the random shift is controlled by
!argument SeparationParam, so that a_new lies withing 
!interval [a_old(1-SeparationParam),a_old(1+SeparationParam)]).

!Arguments:
real(dprec),intent(in)   :: LinDepThreshold
real(dprec),intent(in)   :: SeparationParam
character(Glob_FileNameLength),intent(in) :: FileName
integer,intent(in)       :: PrintInfoSpec

!Local variables:
integer        i,j,k
integer        np,npt,cbs
integer        OpenFileErr,ErrorCode,IFAIL
logical        IsSwapFileOK
integer        BlockSizeForZHEGVX
integer        NumOfEigvalsFound
real(dprec)    Evalue, EVs(1), r
complex(dprec) MaxOverlap,MinOverlap
real(dprec)    AverOverlap
complex(dprec) Min_c,Max_c
real(dprec)    Average_c
integer,allocatable,dimension(:)         :: MaskArray
character(Glob_FileNameLength)          :: ch_temp



if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a39)') 'Routine SeparateLinDepFunc1 has started'
  write(*,'(1x,a27,d23.16)') 'Overlap threshold value is ',LinDepThreshold
  write(*,'(1x,a27,d23.16)') 'Separation paramameter is  ',SeparationParam 
endif 

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
cbs=Glob_CurrBasisSize

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs)) 

!Allocate local workspace
allocate(MaskArray(1:cbs))

!Reading data from swap file

call ReadSwapFileAndDistributeData(IsSwapFileOK)

if (.not.IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
  call ComputeMatElem(1,cbs)
  if (Glob_ProcID==0) write(*,*) ' done'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'


if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in SeparateLinDepFunc1: initial energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*) 'Pair linear dependency check:'
  write(*,*)
endif
do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

!Check overlap
MaskArray(1:cbs)=0
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
    if (abs(Glob_S(j,i))>LinDepThreshold) then
       MaskArray(j)=MaskArray(j)+1
	   MaskArray(i)=MaskArray(j)+1
	   k=k+1
	   if (Glob_ProcID==0) then
	     if (PrintInfoSpec>0) write(*,'(i5,a1,i5,i5,a6,d16.9,a1,d16.9,a6,f11.7)') &
		       k,':',i,j,'   S=(',real(Glob_S(i,j)),',',imag(Glob_S(i,j)),') |S|=',abs(Glob_S(i,j))
         if (PrintInfoSpec==2) then
           write(*,'(i6,1x,<Glob_npt>(1x,d23.16))') i,Glob_NonlinParam(1:Glob_npt,i)
		   write(*,'(6x,a3,d16.9,a1,d16.9,a9,f16.9)') &
		     'c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')    |c|=',abs(Glob_c(i))
           write(*,'(i6,1x,<Glob_npt>(1x,d23.16))') j,Glob_NonlinParam(1:Glob_npt,j)
		   write(*,'(6x,a3,d16.9,a1,d16.9,a9,f16.9)')  &
		     'c=(',real(Glob_c(j)),',',imag(Glob_c(j)),')    |c|=',abs(Glob_c(j))
		 endif
         write(*,*)
       endif
	endif
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients:
Min_c=cmplx(sqrt(huge(ONE)/10),sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) '========== Summary before separation: =========='
endif

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
endif

if (Glob_ProcID==0) then
  write(*,*) 'Basis size before separation ',cbs
  write(*,*) 'Energy before separation     ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif

if (k==0) then
  if (Glob_ProcID==0) then
    write(*,*) 'No linearly dependent functions have been found' 
    write(*,*) 'No file have been written. Program will now stop'
  endif
  stop
endif

j=0
i=1
do i=1,cbs
  if (MaskArray(i)>0) then
    do j=1,npt
      call random_number(r)
	  r=(r-ONEHALF)*2*SeparationParam
      Glob_NonlinParam(j,i)=Glob_NonlinParam(j,i)*(1+r)
    enddo
  endif
enddo
call MPI_BCAST(Glob_NonlinParam,cbs*npt,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
call ComputeMatElem(1,cbs)
if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in EliminateLinDepFunc1: energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*)
  write(*,*) '========== Summary after separation: ==========='
endif

!Check overlap
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients:
Min_c=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
  write(*,*) 'Basis size after separation  ',cbs
  write(*,*) 'Energy after separation      ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif

do i=1,cbs
  Glob_History(i)%Energy=ZERO
  Glob_History(i)%CyclesDone=0
  Glob_History(i)%InitFuncAtLastStep=0
  Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
enddo
Glob_History(cbs)%Energy=Evalue

Glob_CurrEnergy=Evalue
Glob_LastEigvalTol=1.0D30
Glob_BestEigvalTol=1.0D30
Glob_WorstEigvalTol=1.0D-30

ch_temp=Glob_DataFileName
Glob_DataFileName=FileName
if (Glob_ProcID==0) call SaveResults(Sort='no')
Glob_DataFileName=ch_temp

!dellocate local workspace
deallocate(MaskArray)

!deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX) 
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!deallocate global arrays
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) then
  i=len_trim(FileName)
  write(*,*) 'New basis has been saved in file ',FileName(1:i)
  write(*,*) 'Program will now stop'
endif

stop

end subroutine SeparateLinDepFunc1



subroutine SeparateFuncLargeCoeff(LCThreshold,SeparationParam,FileName,PrintInfoSpec)
!Subroutine SeparateFuncLargeCoeff changes linear parameters of
!those basis functions whose linear coefficients (more precisely their
!absulute value) exceed LCThreshold. It is important to note that this 
!subroutine uses normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1). The change
!in nonlinear parameters is controlled by argument SeparationParam.
!Basically the nonlinear parameters of bad functions are chosen
!randomly from the interval [a_old*(1-SeparationParam), a_old*(1+SeparationParam)].
!After saving the results the program is terminated.
!Parameter PrintInfoSpec specifies what information should be
!shown:
! PrintInfoSpec=0 : the subroutine does not print any info
!                   about basis functions.
! PrintInfoSpec=1 : the subroutine prints the values of the linear and 
!                   nonlinear parameters of bad functions.
! PrintInfoSpec=2 : in addition to the latter case it also prints the
!                   linear parameters of all basis functions.

!Arguments:
real(dprec),intent(in)   :: LCThreshold
real(dprec),intent(in)   :: SeparationParam
character(Glob_FileNameLength),intent(in) :: FileName
integer,intent(in)       :: PrintInfoSpec

!Local variables:
integer        i,j,k
integer        np,npt,cbs
integer        OpenFileErr,ErrorCode,IFAIL
logical        IsSwapFileOK
integer        BlockSizeForZHEGVX
integer        NumOfEigvalsFound
real(dprec)    Evalue, EVs(1), r
complex(dprec) MaxOverlap,MinOverlap
real(dprec)    AverOverlap
complex(dprec) Min_c,Max_c
real(dprec)    Average_c
character(Glob_FileNameLength)          :: ch_temp



if (Glob_ProcID==0) then
  write(*,*)
  write(*,'(1x,a42)') 'Routine SeparateFuncLargeCoeff has started'
  write(*,'(1x,a38,d23.16)') 'Linear coefficient threshold value is ',LCThreshold
  write(*,'(1x,a38,d23.16)') 'Separation parameter is               ',SeparationParam 
endif 

!Setting the values of some global variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_HSLeadDim=Glob_CurrBasisSize
np=Glob_np
npt=Glob_npt
Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
cbs=Glob_CurrBasisSize

!Allocate some global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs)) 


!Reading data from swap file

call ReadSwapFileAndDistributeData(IsSwapFileOK)

if (.not.IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
  call ComputeMatElem(1,cbs)
  if (Glob_ProcID==0) write(*,*) ' done'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'


if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in SeparateFuncLargeCoeff: initial energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*)
  if (PrintInfoSpec>=2) then 
    write(*,*) 'Linear coefficients of all basis functions before separation:'
	do i=1,cbs
      write(*,'(1x,i5,a5,d16.9,a1,d16.9,a9,f16.9)')  &
		     i,'  c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')    |c|=',abs(Glob_c(i))
	enddo
  endif
  write(*,*)
endif


!Check overlap
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients
k=0
do i=1,cbs
  if (abs(Glob_c(i))>LCThreshold) then
	k=k+1
    if ((Glob_ProcID==0).and.(PrintInfoSpec>=1)) then
	  if (k==1) write(*,*) 'Functions whose linear coefficients exceed threshold'
      write(*,'(1x,i4,a3,i5,a4,d16.9,a1,d16.9,a8,f16.9)')  &
		     k,':  ',i,'  c=(',real(Glob_c(i)),',',imag(Glob_c(i)),')   |c|=',abs(Glob_c(i))
      write(*,'(1x,<Glob_npt>(1x,d23.16))') Glob_NonlinParam(1:Glob_npt,i)
    endif
  endif
enddo

if (k==0) then 
  if (Glob_ProcID==0) then 
    write(*,*) 'There are no functions whose linear coefficients exceed threshold'
	write(*,*) 'No file have been written. Program will now stop'
  endif
  stop
endif
Min_c=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) '========== Summary before separation: =========='
endif

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
endif

if (Glob_ProcID==0) then
  write(*,*) 'Basis size before separation ',cbs
  write(*,*) 'Energy before separation     ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif

!now change the linear parameters of bad functions
do i=1,cbs
  if (abs(Glob_c(i))>LCThreshold) then
    do j=1,npt
      call random_number(r)
  	  r=(r-ONEHALF)*2*SeparationParam
      Glob_NonlinParam(j,i)=Glob_NonlinParam(j,i)*(1+r)
    enddo
  endif  
enddo
call MPI_BCAST(Glob_NonlinParam,cbs*Glob_npt,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)


if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
call ComputeMatElem(1,cbs)
if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,EVs,Glob_c,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
  Evalue=EVs(1)
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorCode/=0) then
  if (Glob_ProcID==0) write(*,*) &
   'Error in SeparateFuncLargeCoeff: energy cannot be computed'
  stop
endif
Evalue=EVs(1)
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write(*,*) ' done'
  write(*,*)
  write(*,*) '========== Summary after separation: ==========='
endif

!Check overlap
MaxOverlap=cmplx(ZERO,ZERO,dprec)
MinOverlap=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
AverOverlap=ZERO
k=0
do i=1,cbs
  do j=i+1,cbs
	if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
    if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
    AverOverlap=AverOverlap+abs(Glob_S(j,i))
  enddo
enddo
AverOverlap=AverOverlap/(cbs*(cbs-1)/2)

!Check linear coefficients
Min_c=cmplx(sqrt(huge(ONE))/10,sqrt(huge(ONE))/10,dprec)
Max_c=cmplx(ZERO,ZERO,dprec)
Average_c=ZERO
do i=1,cbs
  if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
  if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
  Average_c=Average_c+abs(Glob_c(i))
enddo
Average_c=Average_c/cbs

if (Glob_ProcID==0) then
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal overlap is (',real(MaxOverlap),',',imag(MaxOverlap),')  |S|=',abs(MaxOverlap)
  write(*,'(1x,a20,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal overlap is (',real(MinOverlap),',',imag(MinOverlap),')  |S|=',abs(MinOverlap)
  write(*,'(1x,a37,d16.9)') 'Average absolute value of overlap is ',AverOverlap
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Maximal lin coeff is (',real(Max_c),',',imag(Max_c),')  |c|=',abs(Max_c)
  write(*,'(1x,a22,d16.9,a1,d16.9,a7,f16.9)') &
    'Minimal lin coeff is (',real(Min_c),',',imag(Min_c),')  |c|=',abs(Min_c)
  write(*,'(1x,a41,d16.9)') 'Average absolute value of lin. coeff. is ',Average_c
  write(*,*) 'Basis size after separation  ',cbs
  write(*,*) 'Energy after separation      ',Evalue
  write(*,*) '================================================'
  write(*,*)
endif


do i=1,cbs
  Glob_History(i)%Energy=ZERO
  Glob_History(i)%CyclesDone=0
  Glob_History(i)%InitFuncAtLastStep=0
  Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
enddo
Glob_History(cbs)%Energy=Evalue

Glob_CurrEnergy=Evalue
Glob_LastEigvalTol=1.0D30
Glob_BestEigvalTol=1.0D30
Glob_WorstEigvalTol=1.0D-30

ch_temp=Glob_DataFileName
Glob_DataFileName=FileName
if (Glob_ProcID==0) call SaveResults(Sort='no')
Glob_DataFileName=ch_temp

!deallocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX) 
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!deallocate global arrays
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if (Glob_ProcID==0) then
  i=len_trim(FileName)
  write(*,*) 'New basis has been saved in file ',FileName(1:i)
  write(*,*) 'Program will now stop'
endif

stop

end subroutine SeparateFuncLargeCoeff



subroutine ExpectationValuesG(SymmAdaptMethod,FileName1,FileName2,FileName3,FileName4)
!ExpectationValuesG computes expectation values in the basis of 
!Glob_CurrBasisSize functions. Subroutine DSYGVX is used to solve GSEP.
!Input parameters:
!  SymmAdaptMethod  - defines how the expectation values should be
!calculated. If SymmAdaptMethod=1 then Y^{\dagger}Y operator is applied
!to the ket. This option is faster but requires manual symmetrization
!of the obtained expectation values of such operators as r_{ij} in the
!case of some Hamiltonian symmetries (like charge reversal symmetry). 
!If SymmAdaptMethod=2 then Y operator is applied to both bra and ket. 
!Since the square of the number of terms in Y operator may be significantly 
!larger than the number of terms in Y^{\dagger}Y operator this option is slower. 
!But in this case all the expectation values are correct and ready to use 
!immediately.
!  FileName1 - the name of the file that defines grid for correlation functions.
!              If FileName1 is equal to 'none','NONE', or 'None' then correlation
!              functions are not computed.
!  FileName2 - the name of the file where correlation functions will be stored
!  FileName3 - the name of the file that defines grid for particle densities.
!              If FileName3 is equal to 'none','NONE', or 'None' then particle
!              densities are not computed.
!  FileName4 - the name of the file where particle densities will be stored


!Parameters:
integer,intent(in)  ::    SymmAdaptMethod
character(Glob_FileNameLength),intent(in) :: FileName1
character(Glob_FileNameLength),intent(in) :: FileName2
character(Glob_FileNameLength),intent(in) :: FileName3
character(Glob_FileNameLength),intent(in) :: FileName4

!Local variables:
integer        i,j,k,kk,counter,a,b,c,d
integer        n,np,npt,cbs
integer        OpenFileErr,ErrorCode
logical        IsSwapFileOK
integer        BlockSizeForZHEGVX
integer        NumOfEigvecs,NumOfEigvalsFound
real(dprec)    Evalue
real(dprec),allocatable,dimension(:)      :: Eigvals
complex(dprec),allocatable,dimension(:,:) :: Eigvecs
integer,allocatable,dimension(:)          :: IFAIL
integer        NumOfExpcVals
real(dprec)    factor
real(dprec)    temp1,temp2
integer,allocatable,dimension(:,:)  ::  IdentityPerm
real(dprec)  ReCi,ImCi,ReCj,ImCj,ReOij,ImOij,beta,mu
logical        AreCorrFuncNeeded,ArePartDensNeeded
logical        IsFile1OK,IsFile3OK

!Local variables used to store temporary data
!associated with certain expectation values
integer                                     :: NumCFGridPoints
real(dprec),allocatable,dimension(:)        :: CFGrid
complex(dprec),allocatable,dimension(:,:)   :: CFkl
real(dprec),allocatable,dimension(:,:)      :: CF
integer                                     :: NumDensGridPoints
real(dprec),allocatable,dimension(:)        :: DensGrid
complex(dprec),allocatable,dimension(:,:)   :: Denskl
real(dprec),allocatable,dimension(:,:)      :: Dens
integer                                     :: NumOfCFAndDensExpVals
real(dprec),allocatable,dimension(:)        :: CFDMEkl_s
complex(dprec),allocatable,dimension(:)     :: MEkl
real(dprec),allocatable,dimension(:)        :: MEkl_s, MEkl_e
complex(dprec)                              :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
real(dprec)                                 :: H,S,T,V,MV,drach_MV,Darwin,drach_Darwin,OO
real(dprec)                                 :: H_e,S_e,T_e,V_e,MV_e,drach_MV_e,Darwin_e,drach_Darwin_e,OO_e
complex(dprec),allocatable,dimension(:,:)   :: rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl 
real(dprec),allocatable,dimension(:,:)      :: rm2, rm, r, r2, deltar, drach_deltar
real(dprec),allocatable,dimension(:,:)      :: rm2_e, rm_e, r_e, r2_e, deltar_e, drach_deltar_e 


if (Glob_ProcID==0) then
  write(*,*)
  write(*,*) 'Routine ExpectationValuesG has started'
  write(*,*) 'Number of basis functions',Glob_CurrBasisSize
endif 

!Setting the values of some global and local variables
Glob_GSEPSolutionMethod='G'
Glob_OverlapPenaltyAllowed=.false.
Glob_HSLeadDim=Glob_CurrBasisSize
n=Glob_n
np=Glob_np
npt=Glob_npt
Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
cbs=Glob_CurrBasisSize
NumOfEigvecs=min(cbs,Glob_WhichEigenvalue+10)

!Setting logical variables that determine whether correlation 
!functions and particle densities need be computed
if ((FileName1==' ').or.(FileName1=='none').or. &
    (FileName1=='NONE').or.(FileName1=='None')) then
  AreCorrFuncNeeded=.false.
else
  AreCorrFuncNeeded=.true.  
endif 
if ((FileName3==' ').or.(FileName3=='none').or. &
    (FileName3=='NONE').or.(FileName3=='None')) then
  ArePartDensNeeded=.false.
else
  ArePartDensNeeded=.true.  
endif 

!If corellation functions and/or particle densities are needed
!then we open files FileName1 and FileName3 that contain grids for
!correlation functions and particle densities respectively

!Here we determine the number of grid points for correlation 
!function calculation
if (AreCorrFuncNeeded) then
  if (Glob_ProcID==0) then
    IsFile1OK=.true.
    open(1,file=FileName1,status='old',iostat=OpenFileErr)
    if (OpenFileErr==0) then
      NumCFGridPoints=0
      do while (OpenFileErr==0)
        read (1,*,iostat=OpenFileErr) temp1
        if (temp1<ZERO) then
          OpenFileErr=111
        else  
          NumCFGridPoints=NumCFGridPoints+1
        endif  
      enddo
      NumCFGridPoints=NumCFGridPoints-1
      if ((OpenFileErr>0).or.(NumCFGridPoints==0)) then
        !Improper data in the input file
        write(*,*) 'Error in ExpectationValuesG: improper data in file ',FileName1
        IsFile1OK=.false.
      endif
      close(1)
    else
      IsFile1OK=.false.
    endif
  endif    
  call MPI_BCAST(IsFile1OK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(NumCFGridPoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  !stop if there were problems with correlation function grid file 
  if (.not.IsFile1OK) stop
endif

!Here we determine the number of grid points for particle density
!calculations
if (ArePartDensNeeded) then
  if (Glob_ProcID==0) then
    IsFile3OK=.true.
    open(1,file=FileName3,status='old',iostat=OpenFileErr)
    if (OpenFileErr==0) then
      NumDensGridPoints=0
      do while (OpenFileErr==0)
        read (1,*,iostat=OpenFileErr) temp1
        if (temp1<ZERO) then
          OpenFileErr=111
        else  
          NumDensGridPoints=NumDensGridPoints+1
        endif  
      enddo
      NumDensGridPoints=NumDensGridPoints-1
      if ((OpenFileErr>0).or.(NumDensGridPoints==0)) then
        !Improper data in the input file
        write(*,*) 'Error in ExpectationValuesG: improper data in file ',FileName3
        IsFile3OK=.false.
      endif
      close(1)
    else
      IsFile3OK=.false.
    endif
  endif    
  call MPI_BCAST(IsFile3OK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_BCAST(NumDensGridPoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  !stop if there were problems with particle density grid file 
  if (.not.IsFile3OK) stop
endif

!Allocate arrays that will be used to store grid points and the 
!function values for correlation function and particle density 
!calculations

NumOfCFAndDensExpVals=0

if (AreCorrFuncNeeded) then
  allocate(CFGrid(NumCFGridPoints))
  NumOfCFAndDensExpVals=NumOfCFAndDensExpVals+NumCFGridPoints*n*(n+1)/2
  allocate(CFkl(n*(n+1)/2,NumCFGridPoints))
  allocate(CF(n*(n+1)/2,NumCFGridPoints))
else
  !allocate just one element to have a valid pointer
  allocate(CFGrid(1))
  allocate(CFkl(1,1))  
endif

if (ArePartDensNeeded) then
  allocate(DensGrid(NumDensGridPoints))
  NumOfCFAndDensExpVals=NumOfCFAndDensExpVals+NumDensGridPoints*(n+1)
  allocate(Denskl(n+1,NumDensGridPoints))
  allocate(Dens(n+1,NumDensGridPoints))
else
  !allocate just one element to have a valid pointer
  allocate(DensGrid(1))
  allocate(Denskl(1,1)) 
endif

if((AreCorrFuncNeeded).or.(ArePartDensNeeded)) then
  allocate(CFDMEkl_s(NumOfCFAndDensExpVals))
endif


!Now we open files FileName1 and FileName3 again, but this time
!we read the data from them into arrays CFGrid and DensGrid
if (AreCorrFuncNeeded) then
  if (Glob_ProcID==0) then
    open(1,file=FileName1,status='old')
      do i=1,NumCFGridPoints
        read (1,*) CFGrid(i)
      enddo
    close(1)
  endif  
  call MPI_BCAST(CFGrid,NumCFGridPoints,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
if (ArePartDensNeeded) then
  if (Glob_ProcID==0) then  
    open(1,file=FileName3,status='old')
      do i=1,NumDensGridPoints
        read (1,*) DensGrid(i)
      enddo
    close(1)
  endif  
  call MPI_BCAST(DensGrid,NumDensGridPoints,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

!Allocate global arrays
allocate(Glob_H(cbs,cbs))
allocate(Glob_S(cbs,cbs))
allocate(Glob_diagH(cbs))
allocate(Glob_diagS(cbs))
allocate(Glob_c(cbs))
allocate(Glob_HklBuff1(Glob_HSBuffLen))
allocate(Glob_HklBuff2(Glob_HSBuffLen))
allocate(Glob_SklBuff1(Glob_HSBuffLen))
allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for ZHEGVX
BlockSizeForZHEGVX=ILAENV(1,'ZHETRD','VIU',cbs,cbs,cbs,cbs)
Glob_LWorkForZHEGVX=(BlockSizeForZHEGVX+1)*cbs
allocate(Glob_WorkForZHEGVX(Glob_LWorkForZHEGVX))
allocate(Glob_RWorkForZHEGVX(7*cbs))
allocate(Glob_IWorkForZHEGVX(5*cbs))

!Allocate local arrays

allocate(Eigvals(NumOfEigvecs)) 
allocate(Eigvecs(cbs,NumOfEigvecs))
allocate(IFAIL(cbs))

allocate(IdentityPerm(n,n))
IdentityPerm(1:n,1:n)=0
do i=1,n
  IdentityPerm(i,i)=1
enddo


NumOfExpcVals=6*n*(n+1)/2+9
!rm2          ME(1:n*(n+1)/2)
!rm           ME(n*(n+1)/2+1:2*n*(n+1)/2)
!r            ME(2*n*(n+1)/2+1:3*n*(n+1)/2)
!r2           ME(3*n*(n+1)/2+1:4*n*(n+1)/2)
!deltar       ME(4*n*(n+1)/2+1:5*n*(n+1)/2)
!drach_deltar ME(5*n*(n+1)/2+1:6*n*(n+1)/2)
!H            ME(6*n*(n+1)/2+1)   
!S            ME(6*n*(n+1)/2+2) 
!T            ME(6*n*(n+1)/2+3) 
!V            ME(6*n*(n+1)/2+4) 
!MV           ME(6*n*(n+1)/2+5)
!drach_MVkl   ME(6*n*(n+1)/2+6)
!Darwin       ME(6*n*(n+1)/2+7)
!drach_Darwin ME(6*n*(n+1)/2+8)
!OO           ME(6*n*(n+1)/2+9)   

allocate(MEkl(NumOfExpcVals))
allocate(MEkl_s(NumOfExpcVals))
allocate(MEkl_e(NumOfExpcVals))

allocate(rm2kl(n,n))
allocate(rm2(n,n))
allocate(rm2_e(n,n))

allocate(rmkl(n,n))
allocate(rm(n,n))
allocate(rm_e(n,n))

allocate(rkl(n,n))
allocate(r(n,n))
allocate(r_e(n,n))

allocate(r2kl(n,n))
allocate(r2(n,n))
allocate(r2_e(n,n))

allocate(deltarkl(n,n))
allocate(deltar(n,n))
allocate(deltar_e(n,n))

allocate(drach_deltarkl(n,n))
allocate(drach_deltar(n,n))
allocate(drach_deltar_e(n,n))

call ReadSwapFileAndDistributeData(IsSwapFileOK)

if (.not.IsSwapFileOK) then
  if (Glob_ProcID==0) write(*,'(1x,a52)',advance='no') 'Computing Hamiltonian and overlap matrix elements...'
  call ComputeMatElem(1,cbs)
  if (Glob_ProcID==0) write(*,*) 'done'
endif

do i=1,cbs
  do j=1,i-1
	Glob_H(j,i)=conjg(Glob_H(i,j))
  enddo
  Glob_H(i,i)=Glob_diagH(i)
enddo
do i=1,cbs
  do j=1,i-1
	Glob_S(j,i)=conjg(Glob_S(i,j))
  enddo
  Glob_S(i,i)=cmplx(ONE,ZERO,dprec)
enddo

if (Glob_ProcID==0) then
  write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
  call ZHEGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
     ZERO,ZERO,1,NumOfEigvecs,Glob_AbsTolForZHEGVX, &
	 NumOfEigvalsFound,Eigvals,Eigvecs,cbs,Glob_WorkForZHEGVX,Glob_LWorkForZHEGVX, &
     Glob_RWorkForZHEGVX,Glob_IWorkForZHEGVX,IFAIL,ErrorCode)                 	
  !SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
  !$                   VL, VU, IL, IU, ABSTOL, 
  !$                   M, W, Z, LDZ, WORK, LWORK, 
  !$                   RWORK, IWORK, IFAIL, INFO ) 
endif
call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
if (ErrorCode/=0) then
  if (Glob_ProcID==0) then
    write(*,*) 'failed'
	write(*,*) &
   'Error in ExpectationValuesG: routine ZHEGVX failed with error code ',ErrorCode
  endif
  stop
endif

!sending the eigenvalue and the eigenvector to all processes 
if (Glob_ProcID==0) then
  Evalue=Eigvals(Glob_WhichEigenvalue)
  Glob_c(1:cbs)=Eigvecs(1:cbs,Glob_WhichEigenvalue)
endif
call MPI_BCAST(Evalue,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c,cbs,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
Glob_CurrEnergy=Evalue

!print the lower part of the spectrum
if (Glob_ProcID==0) then
  write(*,*) 'done'
  write(*,*) 'Energy: ',Evalue
  write(*,*)
  write(*,*) 'Lowest eigenvalues:'
  do i=1,NumOfEigvalsFound
    write(*,'(i5,3x,d23.16)') i,Eigvals(i)
  enddo
  write(*,*)
endif

if (Glob_ProcID==0) write(*,'(1x,a31)',advance='no') 'Computing expectation values...'


!main loop
MEkl_s(1:NumOfExpcVals)=ZERO
MEkl_e(1:NumOfExpcVals)=ZERO
CFDMEkl_s(1:NumOfCFAndDensExpVals)=ZERO
counter=0
do i=1,cbs
  do j=1,i  !j=1,i
    counter=counter+1
    if (mod(counter,Glob_NumOfProcs)==Glob_ProcID) then
	  ReCi=real(Glob_c(i),dprec)
	  ImCi=imag(Glob_c(i))
	  ReCj=real(Glob_c(j),dprec)
	  ImCj=imag(Glob_c(j))
      if (i==j) then
	    factor=1/sqrt(Glob_diagS(i)*Glob_diagS(i))  !1
      else
        factor=2/sqrt(Glob_diagS(i)*Glob_diagS(j))  !2
	  endif
	  if (SymmAdaptMethod==1) then
	    do k=1,Glob_NumYHYTerms			    
		  call MatrixElementsForExpcVals(Glob_NonlinParam(1:npt,i),Glob_NonlinParam(1:npt,j), &
		        IdentityPerm,Glob_YHYMatr(1:n,1:n,k),Hkl,Skl,Tkl,Vkl,rm2kl,rmkl,rkl,r2kl,deltarkl,  &
		        drach_deltarkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl,                  &
		        NumCFGridPoints,CFGrid,CFkl,NumDensGridPoints,DensGrid,Denskl,                &
		        AreCorrFuncNeeded,ArePartDensNeeded)	        
		  c=0
          do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=rm2kl(b,a)
            enddo
		  enddo		  
          do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=rmkl(b,a)
            enddo
		  enddo
          do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=rkl(b,a)
            enddo
		  enddo
          do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=r2kl(b,a)
            enddo
		  enddo
          do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=deltarkl(b,a)
            enddo
		  enddo
		  do a=1,n
		    do b=a,n
		      c=c+1; MEkl(c)=drach_deltarkl(b,a)
            enddo
		  enddo
          c=c+1; MEkl(c)=Hkl
          c=c+1; MEkl(c)=Skl  
          c=c+1; MEkl(c)=Tkl
          c=c+1; MEkl(c)=Vkl	
          c=c+1; MEkl(c)=MVkl	
          c=c+1; MEkl(c)=drach_MVkl
          c=c+1; MEkl(c)=Darwinkl	
          c=c+1; MEkl(c)=drach_Darwinkl          
          c=c+1; MEkl(c)=OOkl	                              	  				      

          do a=1,NumOfExpcVals
		    ReOij=real(MEkl(a),dprec)
		    ImOij=imag(MEkl(a))
            beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		    mu=abs(ReCi*(ReCj*ReOij-ImCj*ImOij))+abs(ReCj*(ReCi*ReOij+ImCi*ImOij)) &
		      +abs(ImCi*(ImCj*ReOij+ReCj*ImOij))+abs(ImCj*(ImCi*ReOij-ReCi*ImOij))
            MEkl_s(a)=MEkl_s(a)+factor*Glob_YHYCoeff(k)*beta
		    MEkl_e(a)=MEkl_e(a)+abs(factor*Glob_YHYCoeff(k)*(mu+2*abs(beta)))
		  enddo
		  
		  c=0
		  if (AreCorrFuncNeeded) then
		    do a=1,NumCFGridPoints
		      do b=1,n*(n+1)/2
		        c=c+1 
		        ReOij=real(CFkl(b,a),dprec)
		        ImOij=imag(CFkl(b,a))
		        beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		        CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YHYCoeff(k)*beta 		          
		      enddo
		    enddo
		  endif
		  if (ArePartDensNeeded) then
		    do a=1,NumDensGridPoints
		      do b=1,n+1
		        c=c+1 
		        ReOij=real(Denskl(b,a),dprec)
		        ImOij=imag(Denskl(b,a))
		        beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		        CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YHYCoeff(k)*beta 			          
		      enddo
		    enddo
		  endif 

        enddo !k=1,Glob_NumYHYTerms
      endif !SymmAdaptMethod==1
      

	  if (SymmAdaptMethod==2) then

	    do k=1,Glob_NumYTerms
          do kk=1,Glob_NumYTerms
		    call MatrixElementsForExpcVals(Glob_NonlinParam(1:npt,i),Glob_NonlinParam(1:npt,j),    &
		      Glob_YMatr(1:n,1:n,k),Glob_YMatr(1:n,1:n,kk),Hkl,Skl,Tkl,Vkl,rm2kl,rmkl,rkl,r2kl,deltarkl, &
		      drach_deltarkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl,                          &
              NumCFGridPoints,CFGrid,CFkl,NumDensGridPoints,DensGrid,Denskl,                       &
              AreCorrFuncNeeded,ArePartDensNeeded)
		    c=0
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=rm2kl(b,a)
              enddo
		    enddo		    
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=rmkl(b,a)
              enddo
		    enddo
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=rkl(b,a)
              enddo
		    enddo
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=r2kl(b,a)
              enddo
		    enddo
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=deltarkl(b,a)
              enddo
		    enddo
            do a=1,n
		      do b=a,n
		        c=c+1; MEkl(c)=drach_deltarkl(b,a)
              enddo
		    enddo		    
            c=c+1; MEkl(c)=Hkl
            c=c+1; MEkl(c)=Skl  
            c=c+1; MEkl(c)=Tkl
            c=c+1; MEkl(c)=Vkl	
            c=c+1; MEkl(c)=MVkl	
            c=c+1; MEkl(c)=drach_MVkl
            c=c+1; MEkl(c)=Darwinkl	
            c=c+1; MEkl(c)=drach_Darwinkl          
            c=c+1; MEkl(c)=OOkl	          			      

            do a=1,NumOfExpcVals
		      ReOij=real(MEkl(a),dprec)
		      ImOij=imag(MEkl(a))
              beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		      mu=abs(ReCi*(ReCj*ReOij-ImCj*ImOij))+abs(ReCj*(ReCi*ReOij+ImCi*ImOij)) &
		       +abs(ImCi*(ImCj*ReOij+ReCj*ImOij))+abs(ImCj*(ImCi*ReOij-ReCi*ImOij))
              MEkl_s(a)=MEkl_s(a)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*beta
		      MEkl_e(a)=MEkl_e(a)+abs(factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*(mu+3*abs(beta)))
		    enddo
		    
		    c=0
		    if (AreCorrFuncNeeded) then
		      do a=1,NumCFGridPoints
		        do b=1,n*(n+1)/2
		          c=c+1 
		          ReOij=real(CFkl(b,a),dprec)
		          ImOij=imag(CFkl(b,a))
		          beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		          CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*beta 		          
		        enddo
		      enddo
		    endif
		    if (ArePartDensNeeded) then
		      do a=1,NumDensGridPoints
		        do b=1,n+1
		          c=c+1 
		          ReOij=real(Denskl(b,a),dprec)
		          ImOij=imag(Denskl(b,a))
		          beta=(ReCi*ReCj+ImCi*ImCj)*ReOij+(ImCi*ReCj-ReCi*ImCj)*ImOij
		          CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*beta 			          
		        enddo
		      enddo
		    endif

          enddo !kk=1,Glob_NumYTerms
        enddo !k=1,Glob_NumYTerms
      endif !SymmAdaptMethod==2

    endif
  enddo
enddo
MEkl_e(1:NumOfExpcVals)=MEkl_e(1:NumOfExpcVals)*epsilon(ONE)

!Combining the results of all processes
do a=1,NumOfExpcVals
  temp1=MEkl_s(a)
  call MPI_ALLREDUCE(temp1,temp2,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  MEkl_s(a)=temp2
  temp1=MEkl_e(a)
  call MPI_ALLREDUCE(temp1,temp2,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  MEkl_e(a)=temp2
enddo
k=0
if (AreCorrFuncNeeded) then
  k=NumCFGridPoints*n*(n+1)/2
  call MPI_ALLREDUCE(CFDMEkl_s,CF,k,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
if (ArePartDensNeeded) then
  kk=NumDensGridPoints*(n+1)
  call MPI_ALLREDUCE(CFDMEkl_s(k+1:k+kk),Dens,kk,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)  
endif


!Extracting expectation values from arrays MEkl_s and MEkl_e
c=0
do a=1,n
  do b=a,n
	c=c+1
	rm2(b,a)=MEkl_s(c); rm2(a,b)=MEkl_s(c)
	rm2_e(b,a)=MEkl_e(c); rm2_e(a,b)=MEkl_e(c)
  enddo
enddo
do a=1,n
  do b=a,n
	c=c+1
	rm(b,a)=MEkl_s(c); rm(a,b)=MEkl_s(c)
	rm_e(b,a)=MEkl_e(c); rm_e(a,b)=MEkl_e(c)
  enddo
enddo
do a=1,n
  do b=a,n
	c=c+1
	r(b,a)=MEkl_s(c); r(a,b)=MEkl_s(c)
	r_e(b,a)=MEkl_e(c); r_e(a,b)=MEkl_e(c)
  enddo
enddo
do a=1,n
  do b=a,n
	c=c+1
	r2(b,a)=MEkl_s(c); r2(a,b)=MEkl_s(c)
	r2_e(b,a)=MEkl_e(c); r2_e(a,b)=MEkl_e(c)
  enddo
enddo
do a=1,n
  do b=a,n
	c=c+1
	deltar(b,a)=MEkl_s(c); deltar(a,b)=MEkl_s(c)
	deltar_e(b,a)=MEkl_e(c); deltar_e(a,b)=MEkl_e(c)
  enddo
enddo
do a=1,n
  do b=a,n
	c=c+1
	drach_deltar(b,a)=MEkl_s(c); drach_deltar(a,b)=MEkl_s(c)
	drach_deltar_e(b,a)=MEkl_e(c); drach_deltar_e(a,b)=MEkl_e(c)
  enddo
enddo
c=c+1; H=MEkl_s(c); H_e=MEkl_e(c)
c=c+1; S=MEkl_s(c); S_e=MEkl_e(c)
c=c+1; T=MEkl_s(c); T_e=MEkl_e(c)
c=c+1; V=MEkl_s(c); V_e=MEkl_e(c)
c=c+1; MV=MEkl_s(c); MV_e=MEkl_e(c)
c=c+1; drach_MV=MEkl_s(c); drach_MV_e=MEkl_e(c)
c=c+1; Darwin=MEkl_s(c); Darwin_e=MEkl_e(c)
c=c+1; drach_Darwin=MEkl_s(c); drach_Darwin_e=MEkl_e(c)
c=c+1; OO=MEkl_s(c); OO_e=MEkl_e(c)

!Printing results
if (Glob_ProcID==0) then
  write(*,*) 'done'
  write(*,*) 
  write(*,*) 'Expectation values and the corresponding estimates of'
  write(*,*) 'relative errors due to possible numerical instabilities:'
  write(*,*)
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                     H=',H,'  (',abs(H_e/H),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                     S=',S,'  (',abs(S_e/S),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                     T=',T,'  (',abs(T_e/T),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                     V=',V,'  (',abs(V_e/V),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                    MV=',MV,'  (',abs(MV_e/MV),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '              drach_MV=',drach_MV,'  (',abs(drach_MV_e/drach_MV),')'  
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                  8*MV=',8*MV,'  (',abs(MV_e/MV),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '            8*drach_MV=',8*drach_MV,'  (',abs(drach_MV_e/drach_MV),')'   
  !write(*,'(2x,a23,d23.16)')             '           8*MV(Exact)=',108.1761344D0  
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                Darwin=',Darwin,'  (',abs(Darwin_e/Darwin),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '          drach_Darwin=',drach_Darwin,'  (',abs(drach_Darwin_e/drach_Darwin),')'  
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '                    OO=',OO,'  (',abs(OO_e/OO),')'   
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '          (alpha^2)*MV=',MV*(Glob_FineStructConst**2),'  (',abs(MV_e/MV),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '    (alpha^2)*drach_MV=',drach_MV*(Glob_FineStructConst**2),'  (',abs(drach_MV_e/drach_MV),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '      (alpha^2)*Darwin=',Darwin*(Glob_FineStructConst**2),'  (',abs(Darwin_e/Darwin),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '(alpha^2)*drach_Darwin=',drach_Darwin*(Glob_FineStructConst**2),'  (',abs(drach_Darwin_e/drach_Darwin),')'
  write(*,'(2x,a23,d23.16,a3,d12.5,a1)') '          (alpha^2)*OO=',OO*(Glob_FineStructConst**2),'  (',abs(OO_e/OO),')'      
  write(*,*)
  if ((Glob_NumOfIdentPartSets/=Glob_n+1).and.(SymmAdaptMethod==1)) then
    write(*,*) '(Warning! These values do not account for indistinguishability of'
    write(*,*) 'identical particles and other possible symmetries of the system)'
    write(*,*)
  endif   
  do i=1,n
    write(*,'(11x,a7,i1,a1,d23.16,a3,d12.5,a1)') ' 1/r^2_',i,'=',rm(i,i),'  (',abs(rm_e(i,i)/rm(i,i)),')'
    do j=i+1,n
      write(*,'(11x,a6,i1,i1,a1,d23.16,a3,d12.5,a1)') '1/r^2_',i,j,'=',rm(i,j),'  (',abs(rm_e(i,j)/rm(i,j)),')'
	enddo
  enddo  
  do i=1,n
    write(*,'(11x,a7,i1,a1,d23.16,a3,d12.5,a1)') '   1/r_',i,'=',rm(i,i),'  (',abs(rm_e(i,i)/rm(i,i)),')'
    do j=i+1,n
      write(*,'(11x,a6,i1,i1,a1,d23.16,a3,d12.5,a1)') '  1/r_',i,j,'=',rm(i,j),'  (',abs(rm_e(i,j)/rm(i,j)),')'
	enddo
  enddo
  write(*,*)
  do i=1,n
    write(*,'(11x,a7,i1,a1,d23.16,a3,d12.5,a1)') '     r_',i,'=',r(i,i),'  (',abs(r_e(i,i)/r(i,i)),')'
    do j=i+1,n
      write(*,'(11x,a6,i1,i1,a1,d23.16,a3,d12.5,a1)') '    r_',i,j,'=',r(i,j),'  (',abs(r_e(i,j)/r(i,j)),')'
	enddo
  enddo
  write(*,*)
  do i=1,n
    write(*,'(11x,a7,i1,a1,d23.16,a3,d12.5,a1)') '   r^2_',i,'=',r2(i,i),'  (',abs(r2_e(i,i)/r2(i,i)),')'
    do j=i+1,n
      write(*,'(11x,a6,i1,i1,a1,d23.16,a3,d12.5,a1)') '  r^2_',i,j,'=',r2(i,j),'  (',abs(r2_e(i,j)/r2(i,j)),')'
	enddo
  enddo
  write(*,*)
  do i=1,n
    write(*,'(7x,a10,i1,a2,d23.16,a3,d12.5,a1)') '  delta(r_',i,')=',deltar(i,i),'  (',abs(deltar_e(i,i)/deltar(i,i)),')'
    do j=i+1,n
      write(*,'(7x,a9,i1,i1,a2,d23.16,a3,d12.5,a1)') ' delta(r_',i,j,')=',deltar(i,j),'  (',abs(deltar_e(i,j)/deltar(i,j)),')'
	enddo
  enddo
  write(*,*)  
  do i=1,n
    write(*,'(1x,a16,i1,a2,d23.16,a3,d12.5,a1)') '  drach_delta(r_',i,')=',drach_deltar(i,i),'  (',abs(drach_deltar_e(i,i)/drach_deltar(i,i)),')'
    do j=i+1,n
      write(*,'(1x,a15,i1,i1,a2,d23.16,a3,d12.5,a1)') ' drach_delta(r_',i,j,')=',drach_deltar(i,j),'  (',abs(drach_deltar_e(i,j)/drach_deltar(i,j)),')'
	enddo
  enddo  
  write(*,*)
  if (Glob_NumOfIdentPartSets/=Glob_n+1) then
    write(*,*) 'Based on the particle mass and charge values it was determined'
	write(*,*) 'that the system has the following sets of identical particles:'
    do i=1,Glob_NumOfIdentPartSets
	  j=Glob_NumOfPartInIdentPartSet(i)
      write(*,'(1x,a3,i2,a13,<j>i3)') 'set',i,':   particles',Glob_IdentPartList(1:j,i)
    enddo
	write(*,*) 
    write(*,*) 'Properly symmetrized expectation values of two-particle quantities'
	write(*,*) 'that account for permutational symmetry of the above mentioned sets' 
	write(*,*) 'of identical particles are:'
	write(*,*) '(Warning! An additional symmetrization might be necessary if the'
	write(*,*) 'Young operator contains other types of symmetries)' 
	write(*,*)
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a6,i1,a3)',advance='no') '1/r^2_',a,' = '
        else
          write(*,'(a6,i1,i1,a3)',advance='no') '1/r^2_',a,b,' = '
		endif
		beta=beta+rm2(a,b)
		mu=mu+rm2_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo	
	write(*,*)	
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a4,i1,a3)',advance='no') '1/r_',a,' = '
        else
          write(*,'(a4,i1,i1,a3)',advance='no') '1/r_',a,b,' = '
		endif
		beta=beta+rm(a,b)
		mu=mu+rm_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo
	write(*,*)
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a2,i1,a3)',advance='no') 'r_',a,' = '
        else
          write(*,'(a2,i1,i1,a3)',advance='no') 'r_',a,b,' = '
		endif
		beta=beta+r(a,b)
		mu=mu+r_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo
	write(*,*)
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a4,i1,a3)',advance='no') 'r^2_',a,' = '
        else
          write(*,'(a4,i1,i1,a3)',advance='no') 'r^2_',a,b,' = '
		endif
		beta=beta+r2(a,b)
		mu=mu+r2_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo
	write(*,*)
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a8,i1,a4)',advance='no') 'delta(r_',a,') = '
        else
          write(*,'(a8,i1,i1,a4)',advance='no') 'delta(r_',a,b,') = '
		endif
		beta=beta+deltar(a,b)
		mu=mu+deltar_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo
	write(*,*)
    do i=1,Glob_NumOfNoneqvPairSets
	  beta=ZERO
	  mu=ZERO
	  k=Glob_NumOfPairsInEqvPairSet(i)
	  write(*,'(1x)',advance='no')
	  do j=1,k
	    a=Glob_EqvPairList(1,j,i)
		b=Glob_EqvPairList(2,j,i)
	    if (a==b) then
          write(*,'(a14,i1,a4)',advance='no') 'drach_delta(r_',a,') = '
        else
          write(*,'(a14,i1,i1,a4)',advance='no') 'drach_delta(r_',a,b,') = '
		endif
		beta=beta+drach_deltar(a,b)
		mu=mu+drach_deltar_e(a,b)
	  enddo
      write(*,'(d23.16,a2,d12.5,a1)') beta/k,' (',abs(mu/(beta)),')'
    enddo
	write(*,*)	
  endif

    !write(*,*)
  	!write(*,*) 'Glob_NumOfNoneqvPairSets=',Glob_NumOfNoneqvPairSets
  	!write(*,*)
  	!do i=1,Glob_NumOfNoneqvPairSets
  	!  write(*,'(1x,a28,i1,a2,i2)') 'Glob_NumOfPairsInEqvPairSet(',i,')=',Glob_NumOfPairsInEqvPairSet(i)
  	!enddo
  	!write(*,*)
  	!do i=1,Glob_NumOfNoneqvPairSets
  	!  k=Glob_NumOfPairsInEqvPairSet(i)
  	!  do j=1,k
  	!    write(*,'(1x,a19,i1,a1,i1,a2,i2)') 'Glob_EqvPairList(1,',j,',',i,')=',Glob_EqvPairList(1,j,i)
  	!    write(*,'(1x,a19,i1,a1,i1,a2,i2)') 'Glob_EqvPairList(2,',j,',',i,')=',Glob_EqvPairList(2,j,i)
  	!  enddo
  	!  write(*,*)
  	!enddo  

  !Storing correlation functions
  if (AreCorrFuncNeeded) then
    open(1,file=FileName2,status='replace')
    !first we print titles of all data columns
    write(1,'(11x,a2)',advance='no') 'r '
    do i=1,Glob_NumOfNoneqvPairSets
      a=Glob_EqvPairList(1,1,i)
      b=Glob_EqvPairList(2,1,i)
	  if (a>b) then
		c=a; a=b; b=c
	  endif      
      if (a==b) then
        write(1,'(21x,a1,i1,1x)',advance='no') 'g',a
      else
        write(1,'(21x,a1,i1,i1)',advance='no') 'g',a,b
      endif
    enddo
    write(1,*)
    !then we print the data columns themselves
    do kk=1,NumCFGridPoints
      write(1,'(1x,d23.16)',advance='no') CFGrid(kk)
      do i=1,Glob_NumOfNoneqvPairSets
	    mu=ZERO
	    k=Glob_NumOfPairsInEqvPairSet(i)
	    do j=1,k
	      a=Glob_EqvPairList(1,j,i)
		  b=Glob_EqvPairList(2,j,i)
		  if (a>b) then
		    c=a; a=b; b=c
		  endif
		  mu=mu+CF(b-(a-1)*(a-2*n)/2,kk)	    
	    enddo
	    write(1,'(1x,d23.16)',advance='no') mu/k 
      enddo
      write(1,*)
    enddo  
    close(1)
    i=len_trim(FileName2)
    write(*,*) 'Correlation functions have been stored in file ',FileName2(1:i)
    write(*,*)
  endif
  
  !Storing particle densisities
  if (ArePartDensNeeded) then
    open(1,file=FileName4,status='replace')
    !first we print titles of all data columns
    write(1,'(11x,a3)',advance='no') 'r  '
    do i=1,Glob_NumOfIdentPartSets
      write(1,'(20x,a3,i1)',advance='no') 'rho',Glob_IdentPartList(1,i)
    enddo   
    write(1,*) 
    !then we print the data columns themselves
    do kk=1,NumDensGridPoints
      write(1,'(1x,d23.16)',advance='no') DensGrid(kk)
      do i=1,Glob_NumOfIdentPartSets
        mu=ZERO
        k=Glob_NumOfPartInIdentPartSet(i)
        do j=1,k
          mu=mu+Dens(Glob_IdentPartList(j,i),kk)
        enddo
        write(1,'(1x,d23.16)',advance='no') mu/k
      enddo
      write(1,*)
    enddo       
    close(1)
    i=len_trim(FileName4)
    write(*,*) 'Particle densities have been stored in file ',FileName4(1:i)
    write(*,*)  
  endif

endif

!deallocate local arrays

deallocate(rm2kl)
deallocate(rm2)
deallocate(rm2_e)

deallocate(rmkl)
deallocate(rm)
deallocate(rm_e)

deallocate(rkl)
deallocate(r)
deallocate(r_e)

deallocate(r2kl)
deallocate(r2)
deallocate(r2_e)

deallocate(deltarkl)
deallocate(deltar)
deallocate(deltar_e)

deallocate(drach_deltarkl)
deallocate(drach_deltar)
deallocate(drach_deltar_e)

deallocate(MEkl)
deallocate(MEkl_s)
deallocate(MEkl_e)

deallocate(IdentityPerm)

deallocate(IFAIL)
deallocate(Eigvecs)
deallocate(Eigvals)

!dellocate workspace for ZHEGVX
deallocate(Glob_IWorkForZHEGVX) 
deallocate(Glob_RWorkForZHEGVX)
deallocate(Glob_WorkForZHEGVX)

!deallocate global arrays
deallocate(Glob_SklBuff2)
deallocate(Glob_SklBuff1)
deallocate(Glob_HklBuff2)
deallocate(Glob_HklBuff1)
deallocate(Glob_c)
deallocate(Glob_diagS)
deallocate(Glob_diagH)
deallocate(Glob_S)
deallocate(Glob_H)

if((AreCorrFuncNeeded).or.(ArePartDensNeeded)) then
  deallocate(CFDMEkl_s)
endif

deallocate(DensGrid) 
deallocate(Denskl)  
if (ArePartDensNeeded) deallocate(Dens)  
deallocate(CFGrid)
deallocate(CFkl) 
if (AreCorrFuncNeeded) deallocate(CF)
 

if (Glob_ProcID==0) write (*,*) 'Routine ExpectationValuesG has finished'

end subroutine ExpectationValuesG




end module workproc
