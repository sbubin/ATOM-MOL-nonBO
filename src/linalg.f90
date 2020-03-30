!DEC$ DECLARE

module linalg !version 0.52

use globvars

!Module linalg contains some linear algebra routines
!that are used in calculations

!These global variables define in what mode (parallel or using 
!a single process only) certain routines should work. They
!also define whether BLAS or generic routines should be
!used (it concernes only the single process mode). The
!values of these variables may change depending on the size
!of the eigenvalue problem, architechture, etc.
!Basically, "0" stands for "no" and "1" stands for "yes".
integer,parameter :: UseBlas=1
integer,parameter :: PMode=1
integer,parameter :: SMode=1

integer :: Glob_LDLTF_PMode=0
integer :: Glob_LDLTF_SModeType=0

integer :: Glob_LDLHF_PMode=0
integer :: Glob_LDLHF_SModeType=0

integer :: Glob_LDLTS_PMode=PMode
integer :: Glob_LDLTS_SModeType=SMode
integer :: Glob_LDLTS_UseBLAS=UseBlas

integer :: Glob_LDLHS_PMode=PMode
integer :: Glob_LDLHS_SModeType=SMode
integer :: Glob_LDLHS_UseBLAS=UseBlas

integer :: Glob_MTMVL_PMode=PMode
integer :: Glob_MTMVL_UseBLAS=UseBlas

integer :: Glob_MHMVL_PMode=PMode
integer :: Glob_MHMVL_UseBLAS=UseBlas

integer :: Glob_MTMV_PMode=PMode
integer :: Glob_MTMV_UseBLAS=UseBlas

integer :: Glob_MHMV_PMode=PMode
integer :: Glob_MHMV_UseBLAS=UseBlas

integer :: Glob_VMMTMV_PMode=PMode

integer :: Glob_VMMHMV_PMode=PMode

integer :: Glob_RMaxAbsEl_PMode=PMode

integer :: Glob_CMaxAbsReOrIm_PMode=PMode

integer :: Glob_RDotProd_PMode=PMode
integer :: Glob_RDotProd_UseBLAS=UseBlas

integer :: Glob_CDotProd_PMode=PMode
integer :: Glob_CDotProd_UseBLAS=UseBlas

integer :: Glob_RDotProdItself_PMode=PMode
integer :: Glob_RDotProdItself_UseBLAS=UseBlas

integer :: Glob_CDotProdItself_PMode=PMode
integer :: Glob_CDotProdItself_UseBLAS=UseBlas

integer :: Glob_RDotProdQuotient_PMode=PMode

integer :: Glob_CDotProdQuotient_PMode=PMode

integer :: Glob_RVScale_UseBLAS=UseBlas

integer :: Glob_CVScale_UseBLAS=UseBlas

integer :: Glob_RVDiffEucNorm_PMode=PMode

integer :: Glob_CVDiffEucNorm_PMode=PMode


contains


subroutine linalg_setparam(n)
!Subroutine linalg_setparam sets the values of the global parameters
!declared in module linalg that provide maximum performance (roughly).
!The user needs to pass the dimension of the linear algebra problems, n.
!The optimal values of the global parameters are specific for each 
!particular n. Therefore it is important that the user calls linalg_setparam
!each time the dimensionality of the linear algebra problem is changed
!in his/her calculations.  
!The optimal values of the global parameters also quite depend on the 
!computer architecture and the processor interconnect speed. Below one can
!find an implementation that is based on rough tests made on a SGI Altix 4700 
!supercomputer and on an a AMD Athlon MP PC cluster with Gigabit Ethernet.
!Each time the present code is used on a different machine, it is highly
!recommended to run performance tests and determine the sets of optimal 
!parameters as a function of n.
character(20),parameter :: machine='marin'   !marin,beowulf
integer n

select case (machine)

case('marin') !Optimal parameters for SGI Altix 4700
  !Glob_LDLTF_PMode, Glob_LDLTF_SModeType
  select case (n)
  case(1:1100)
    Glob_LDLTF_PMode=0; Glob_LDLTF_SModeType=1;  
  case(1101:)
    Glob_LDLTF_PMode=0; Glob_LDLTF_SModeType=0;
  endselect
  if ((Glob_NumOfProcs>=8).and.(n>6000)) Glob_LDLTF_PMode=1
  if ((Glob_NumOfProcs>=16).and.(n>5000)) Glob_LDLTF_PMode=1
  !Glob_LDLHF_PMode, Glob_LDLHF_SModeType
  select case (n)
  case(1:500)
    Glob_LDLHF_PMode=0; Glob_LDLHF_SModeType=1;  
  case(501:)
    Glob_LDLHF_PMode=0; Glob_LDLHF_SModeType=0;  
  endselect
  if ((Glob_NumOfProcs>=8).and.(n>5000)) Glob_LDLHF_PMode=1
  if ((Glob_NumOfProcs>=16).and.(n>3000)) Glob_LDLHF_PMode=1  

case('beowulf') !Optimal parameters for AMD Athlon MP cluster with Gigabit Ethernet. 
  !Glob_LDLTF_PMode, Glob_LDLTF_SModeType
  select case (n)
  case(1:1100)
    Glob_LDLTF_PMode=0; Glob_LDLTF_SModeType=1;  
  case(1101:)
    Glob_LDLTF_PMode=0; Glob_LDLTF_SModeType=0;
  endselect
  if ((Glob_NumOfProcs>=8).and.(n>6000)) Glob_LDLTF_PMode=1
  if ((Glob_NumOfProcs>=16).and.(n>5000)) Glob_LDLTF_PMode=1
  !Glob_LDLHF_PMode, Glob_LDLHF_SModeType
  select case (n)
  case(1:500)
    Glob_LDLHF_PMode=0; Glob_LDLHF_SModeType=1;  
  case(501:)
    Glob_LDLHF_PMode=0; Glob_LDLHF_SModeType=0;  
  endselect
  if ((Glob_NumOfProcs>=8).and.(n>5000)) Glob_LDLHF_PMode=1
  if ((Glob_NumOfProcs>=16).and.(n>3000)) Glob_LDLHF_PMode=1  

case default

endselect

end subroutine linalg_setparam


subroutine LDLTF(m,n,A,nA,invD,w,ErrorCode)
!Subroutine LDLTF upadates the factorization of a real symmetric 
!matrix A of size n in the form A=L*D*LT using a known factorization 
!of its submatrix of size m-1. Here L is a lower triangular matrix 
!with diagonal elements equal to 1, D is a diagonal matrix, LT is 
!the transpose of L. If no factorization of a submatrix is 
!known, one should set m=1. The lower triangle of A (including the 
!diagonal) remains unchanged on exit. 
!  Input parameters :
!     m-1 - The size of the submatrix of A, whose factorization is 
!           known (m>=1);
!       n - The size of matrix A;
!      nA - The leading dimension of A;
!       A - A two-dimensional array. Its lower triangle (including 
!           the diagonal) contains matix A and the upper one
!           contains the elements of matrix LT (of size m-1);
!    invD - An array containing inversed values of
!           the diagonal elements of D (up to size m-1);
!       w - A work vector (at least of size n)
!  Output parameters:
!       A - A two-dimensional array. The elements of matrix A 
!           (including the diagonal) are stored in the lower 
!           triangle. Matrix LT is returned in the upper 
!           triangle (updated to size n);
!    invD - An array containing inversed values of
!           the diagonal elements of D (updated to size n);
!ErrorCode- The error flag. If ErrorCode=0 then the procedure
!           finished succesfully. If ErrorCode=1 then the
!           procedure cannot be completed because matrix A is 
!           singular

!Arguments :  
integer	       m,n,nA,ErrorCode
real(dprec)    A(nA,n),invD(n),w(n)
!Local variables :
integer        i,j,jm,im
integer        q,k,ji,jf,jim,jiR,RowsPerProc,mod_im_Glob_NumOfProcs
real(dprec)    x,y,z

ErrorCode=0
if (Glob_LDLTF_PMode==0) then
  !Single process version
  if ((Glob_LDLTF_SModeType/=0).or.(Glob_ProcID==0)) then
    do i=m,n
      do j=1,i      
        jm=j-1 
        if (i==j) then
          !DEC$ UNROLL (4)
          do q=1,jm
            w(q)=A(q,i)*invD(q)
          enddo
	      x=A(i,j)
	      !DEC$ UNROLL (4)
	      do q=1,jm
	        x=x-A(q,i)*w(q)
	      enddo
	      !DEC$ UNROLL (4)      
	      do q=1,jm
	        A(q,i)=w(q)
	      enddo					
	      if (x==ZERO) then 
	  	    ErrorCode=1
		    return
	      endif
	      invD(j)=ONE/x	 
	    else
	      x=A(i,j)
	      !DEC$ UNROLL (4)
	      do q=1,jm
	        x=x-A(q,i)*A(q,j)
	      enddo	
	      A(j,i)=x    
	    endif
      enddo
    enddo
  endif
  if (Glob_LDLTF_SModeType==0) then
    do i=m,n
      call MPI_BCAST(A(1:i-1,i),i-1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    enddo
    call MPI_BCAST(invD(m:n),n-m+1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  endif  
else
  !Parallel version (far from being perfect...)
  do i=m,n
    im=i-1
	A(1:im,i)=ZERO
	RowsPerProc=im/Glob_NumOfProcs
	mod_im_Glob_NumOfProcs=mod(im,Glob_NumOfProcs)
	do k=1,RowsPerProc
      jf=k*Glob_NumOfProcs
	  jim=jf-Glob_NumOfProcs
	  ji=jim+1
	  jiR=ji+Glob_ProcID
	  !A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
	  z=ZERO
	  !DEC$ UNROLL (4)
	  do q=1,jim
	    z=z-A(q,i)*A(q,jiR)
	  enddo
	  A(jiR,i)=z
	  !
      call MPI_ALLREDUCE(A(ji:jf,i),w(1:Glob_NumOfProcs),Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      A(ji:jf,i)=w(1:Glob_NumOfProcs)
	  do j=ji,jf
	    jm=j-1
		!A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
		z=A(j,i)+A(i,j)
		!DEC$ UNROLL (4)
		do q=ji,jm
		  z=z-A(q,i)*A(q,j)
		enddo
		A(j,i)=z
		!
	  enddo
	enddo
    if (mod_im_Glob_NumOfProcs>0) then
      ji=RowsPerProc*Glob_NumOfProcs+1
	  jim=ji-1
	  jf=im
	  jiR=ji+Glob_ProcID
      if (jiR<i) then
	    !A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
	    z=ZERO
	    !DEC$ UNROLL (4)
	    do q=1,jim
	      z=z-A(q,i)*A(q,jiR)
	    enddo
	    A(jiR,i)=z
	    !
	  endif
      call MPI_ALLREDUCE(A(ji:jf,i),w(1:mod_im_Glob_NumOfProcs),mod_im_Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      A(ji:jf,i)=w(1:mod_im_Glob_NumOfProcs)
	  do j=ji,jf
        jm=j-1
		!A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
		z=A(j,i)+A(i,j)
		!DEC$ UNROLL (4)
		do q=ji,jm
		  z=z-A(q,i)*A(q,j)
		enddo
		A(j,i)=z
		!
	  enddo
	endif
	!j==1 case
	w(1:im)=A(1:im,i)*invD(1:im)
    y=ZERO
	do k=1+Glob_ProcID,im,Glob_NumOfProcs
       y=y+A(k,i)*w(k)
	enddo
    call MPI_ALLREDUCE(y,x,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    x=A(i,i)-x
	A(1:im,i)=w(1:im)
	if (x==ZERO) then
      ErrorCode=1
	  return
	endif
	invD(i)=ONE/x
  enddo
endif

end subroutine LDLTF


subroutine LDLHF(m,n,A,nA,invD,w,ErrorCode)
!Subroutine LDLHF upadates the factorization of a complex hermitian 
!matrix A of size n in the form A=L*D*LH using a known factorization 
!of its submatrix of size m-1. Here L is a lower triangular matrix 
!with diagonal elements equal to 1, D is a diagonal matrix, LH is 
!the hermitian conjugate of L. If no factorization of a submatrix is 
!known, one should set m=1. The lower triangle of A (including the 
!diagonal) remains unchanged on exit. 
!  Input parameters :
!     m-1 - The size of the submatrix of A, whose factorization is 
!           known (m>=1);
!       n - The size of matrix A;
!      nA - The leading dimension of A;
!       A - A two-dimensional array. Its lower triangle (including 
!           the diagonal) contains matix A and the upper one
!           contains the elements of matrix LH (of size m-1);
!    invD - An array containing inversed values of
!           the diagonal elements of D (up to size m-1);
!       w - A work vector (at least of size n)
!  Output parameters:
!       A - A two-dimensional array. The elements of matrix A 
!           (including the diagonal) are stored in the lower 
!           triangle. Matrix LH is returned in the upper 
!           triangle (updated to size n);
!    invD - An array containing inversed values of
!           the diagonal elements of D (updated to size n);
!ErrorCode- The error flag. If ErrorCode=0 then the procedure
!           finished succesfully. If ErrorCode=1 then the
!           procedure cannot be completed because matrix A is 
!           singular

!Arguments :  
integer	       m,n,nA,ErrorCode
complex(dprec) A(nA,n),invD(n),w(n)
!Local variables :
integer        i,j,jm,im
integer        q,k,ji,jf,jim,jiR,RowsPerProc,mod_im_Glob_NumOfProcs
complex(dprec) x,y

ErrorCode=0
if (Glob_LDLHF_PMode==0) then
  !Single process version
  if ((Glob_LDLHF_SModeType/=0).or.(Glob_ProcID==0)) then
    do i=m,n
      do j=1,i      
        jm=j-1 
        if (i==j) then
          !DEC$ UNROLL (4)
          do q=1,jm
            w(q)=A(q,i)*invD(q)
          enddo
	      x=A(i,j)
	      !DEC$ UNROLL (4)
	      do q=1,jm
	        x=x-conjg(A(q,i))*w(q)
	      enddo  
	      x=conjg(x)
	      !DEC$ UNROLL (4)
	      do q=1,jm
	        A(q,i)=w(q)
	      enddo					
	      if (x==ZERO) then 
	  	    ErrorCode=1
		    return
	      endif
	      invD(j)=ONE/x	 
	    else
	      x=A(i,j)
	      !DEC$ UNROLL (4)
	      do q=1,jm
	        x=x-conjg(A(q,i))*A(q,j)
	      enddo	    
	      A(j,i)=conjg(x)	
	    endif
      enddo
    enddo
  endif
  if (Glob_LDLHF_SModeType==0) then
    do i=m,n
      call MPI_BCAST(A(1:i-1,i),i-1,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    enddo
    call MPI_BCAST(invD(m:n),n-m+1,MPI_DCOMP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  endif   
else
  !Parallel version (far from being perfect...)
  do i=m,n
    im=i-1
	A(1:im,i)=ZERO
	RowsPerProc=im/Glob_NumOfProcs
	mod_im_Glob_NumOfProcs=mod(im,Glob_NumOfProcs)
	do k=1,RowsPerProc
      jf=k*Glob_NumOfProcs
	  jim=jf-Glob_NumOfProcs
	  ji=jim+1
	  jiR=ji+Glob_ProcID
	  A(jiR,i)=conjg(-dot_product(A(1:jim,i),A(1:jim,jiR)))
      call MPI_ALLREDUCE(A(ji:jf,i),w(1:Glob_NumOfProcs),Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      A(ji:jf,i)=w(1:Glob_NumOfProcs)
      !DEC$ UNROLL (4)
	  do j=ji,jf
	    jm=j-1
		A(j,i)=A(j,i)+conjg(A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j)))
	  enddo
	enddo
    if (mod_im_Glob_NumOfProcs>0) then
      ji=RowsPerProc*Glob_NumOfProcs+1
	  jim=ji-1
	  jf=im
	  jiR=ji+Glob_ProcID
      if (jiR<i) then
	    A(jiR,i)=conjg(-dot_product(A(1:jim,i),A(1:jim,jiR)))
	  endif
      call MPI_ALLREDUCE(A(ji:jf,i),w(1:mod_im_Glob_NumOfProcs),mod_im_Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      A(ji:jf,i)=w(1:mod_im_Glob_NumOfProcs)
      !DEC$ UNROLL (4)
	  do j=ji,jf
        jm=j-1
		A(j,i)=A(j,i)+conjg(A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j)))
	  enddo
	endif
	!j==1 case
	w(1:im)=A(1:im,i)*invD(1:im)
    y=ZERO
	do k=1+Glob_ProcID,im,Glob_NumOfProcs
       y=y+conjg(A(k,i))*w(k)
	enddo
	y=conjg(y)
    call MPI_ALLREDUCE(y,x,1,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    x=A(i,i)-x
	A(1:im,i)=w(1:im)
	if (x==ZERO) then
      ErrorCode=1
	  return
	endif
	invD(i)=ONE/x
  enddo
endif

end subroutine LDLHF


subroutine LDLTS(n,A,nA,invD,b,x)
!Subroutine LDLTS finds the solution of a linear system
!A*x=b, where A is a real symmetric matrix, b is a right-hand 
!side. Prior to calling LDLTS one must factorize A in the
!form A=L*D*LT using routine LDLTF. The solution of the 
!system is found as a result of consecutive solutions of
!systems L*y=b and D*LT*x=y. Vector b is destroyed on exit (only
!if the routine runs in parallel mode). 
!  Input parameters :
!       n - The size of the linear system;
!      nA - the leading dimension of A;
!       A - A two-dimensional array that contains the elements
!           of matrix LT (the result of routine LDLTF) in the upper 
!           triangle (not including the diagonal elements, which 
!           are assumed to be ones). The lower part of A is not used.
!    invD - An array containing inversed values of
!           the diagonal elements of D (the result of LDLTF);
!       b - An array containing the right-hand side. It
!           is also used as workspace so it is destroyed
!           on exit.
!  Output parameters :
!       x - An array containing the solution;

!Arguments
integer     n,nA
real(dprec) A(nA,n),invD(n),b(n),x(n)
!Local variables
integer     i,j,k,ji,jf,jim,jfm,RowsPerProc,mod_n_Glob_NumOfProcs,jiR
real(dprec) t

if (Glob_LDLTS_PMode==0) then
  if (Glob_LDLTS_UseBLAS==0) then
    !Solution of L*y=b
    x(1:n)=b(1:n)
    do j=1,n
      t=x(j)
      !DEC$ UNROLL (4)
      do i=1,j-1
        t=t-A(i,j)*x(i)
      enddo
      x(j)=t
    enddo  
    !Solution of D*LT*x=y  
    x(1:n)=x(1:n)*invD(1:n) 
    do j=N,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
      do i=j-1,1,-1
        x(i)=x(i)-t*A(i,j)
      enddo
    enddo         
  else
    x(1:n)=b(1:n)
    !call BLAS routine DTRSV   
    call DTRSV('U','T','U',n,A,nA,x,1)
    x(1:n)=x(1:n)*invD(1:n)
    !call BLAS routine DTRSV
    call DTRSV('U','N','U',n,A,nA,x,1)
  endif	
else
  !Solution of L*y=b
  RowsPerProc=n/Glob_NumOfProcs
  x(1:n)=ZERO
  do i=1,RowsPerProc
    ji=(i-1)*Glob_NumOfProcs+1
	jim=ji-1
    jf=jim+Glob_NumOfProcs
	t=ZERO
	jiR=ji+Glob_ProcID
	!DEC$ UNROLL (4)
	do k=1,jim
	  t=t-A(k,jiR)*x(k)
	enddo
	x(jiR)=t+b(jiR)
    call MPI_ALLREDUCE(x(ji:jf),b(ji:jf),Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	x(ji:jf)=b(ji:jf)
	do j=ji,jf
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=ji,j-1
	    t=t-A(k,j)*x(k)
      enddo
	  x(j)=t
	enddo
  enddo
  mod_n_Glob_NumOfProcs=mod(n,Glob_NumOfProcs)
  if (mod_n_Glob_NumOfProcs>0) then
    jim=RowsPerProc*Glob_NumOfProcs
	ji=jim+1
	jiR=ji+Glob_ProcID
	if (jiR<=n) then
	  t=ZERO
	  !DEC$ UNROLL (4)
      do k=1,jim
        t=t-A(k,jiR)*x(k)
	  enddo
	  x(jiR)=t+b(jiR)
    endif
    call MPI_ALLREDUCE(x(ji:n),b(ji:n),mod_n_Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    x(ji:n)=b(ji:n)
	do j=ji,n
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=ji,j-1
	    t=t-A(k,j)*x(k)
      enddo
	  x(j)=t
	enddo
  endif
  !Solution of D*LT*x=y  
  b(1:n)=x(1:n)*invD(1:n)   
  x(1:n)=ZERO
  do i=1,RowsPerProc
    jf=n-i*Glob_NumOfProcs+1
	ji=jf-1+Glob_NumOfProcs
	jiR=ji+1+Glob_ProcID
	if (jiR<=n) then
 	  t=x(jiR)
 	  !DEC$ UNROLL (4)
      do k=ji,1,-1
         x(k)=x(k)-t*A(k,jiR)
	  enddo
      call MPI_ALLREDUCE(x(jf:ji),b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs),Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      x(jf:ji)=b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
	endif
    x(jf:ji)=x(jf:ji)+b(jf:ji)
	do j=ji,jf,-1
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=j-1,jf,-1
	     x(k)=x(k)-t*A(k,j)
	  enddo
	enddo
  enddo
  if (mod_n_Glob_NumOfProcs>0) then
	ji=mod_n_Glob_NumOfProcs
	jiR=ji+1+Glob_ProcID
	if (jiR<=n) then
 	  t=x(jiR)
 	  !DEC$ UNROLL (4)
      do k=ji,1,-1
         x(k)=x(k)-t*A(k,jiR)
	  enddo
      call MPI_ALLREDUCE(x(1:ji),b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs),mod_n_Glob_NumOfProcs,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      x(1:ji)=b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
	endif
	x(1:ji)=x(1:ji)+b(1:ji)
	do j=ji,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=j-1,1,-1
	     x(k)=x(k)-t*A(k,j)
	  enddo
	enddo
  endif
endif

end subroutine LDLTS


subroutine LDLHS(n,A,nA,invD,b,x)
!Subroutine LDLHS finds the solution of a linear system
!A*x=b, where A is a complex hermitian matrix, b is a right-hand 
!side. Prior to calling LDLHS one must factorize A in the
!form A=L*D*LH using routine LDLHF. The solution of the 
!system is found as a result of consecutive solutions of
!systems L*y=b and D*LH*x=y. Vector b is destroyed on exit (only
!if the routine runs in parallel mode). 
!  Input parameters :
!       n - The size of the linear system;
!      nA - the leading dimension of A;
!       A - A two-dimensional array that contains the elements
!           of matrix LH (the result of routine LDLHF) in the upper 
!           triangle (not including the diagonal elements, which 
!           are assumed to be ones). The lower part of A is not used.
!    invD - An array containing inversed values of
!           the diagonal elements of D (the result of LDLHF);
!       b - An array containing the right-hand side. It
!           is also used as workspace so it is destroyed
!           on exit.
!  Output parameters :
!       x - An array containing the solution;

!Arguments
integer        n,nA
complex(dprec) A(nA,n),invD(n),b(n),x(n)
!Local variables
integer        i,j,k,ji,jf,jim,jfm,RowsPerProc,mod_n_Glob_NumOfProcs,jiR
complex(dprec) t

if (Glob_LDLHS_PMode==0) then
  if (Glob_LDLHS_UseBLAS==0) then
    !Solution of L*y=b
    x(1:n)=b(1:n)
    do j=1,n
      t=x(j)
      !DEC$ UNROLL (4)
      do i=1,j-1
        t=t-conjg(A(i,j))*x(i)
      enddo
      x(j)=t
    enddo  	
    !Solution of D*LH*x=y  
    x(1:n)=x(1:n)*invD(1:n) 
    do j=N,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
      do i=j-1,1,-1
        x(i)=x(i)-t*A(i,j)
      enddo
    enddo  
  else
    x(1:n)=b(1:n)
    !call BLAS routine ZTRSV
    call ZTRSV('U','C','U',n,A,nA,x,1)
    x(1:n)=x(1:n)*invD(1:n)
    !call BLAS routine ZTRSV
    call ZTRSV('U','N','U',n,A,nA,x,1)
  endif
else
  !Solution of L*y=b
  RowsPerProc=n/Glob_NumOfProcs
  x(1:n)=ZERO
  do i=1,RowsPerProc
    ji=(i-1)*Glob_NumOfProcs+1
	jim=ji-1
    jf=jim+Glob_NumOfProcs
	t=ZERO
	jiR=ji+Glob_ProcID
	!DEC$ UNROLL (4)
	do k=1,jim
	  t=t-conjg(A(k,jiR))*x(k)
	enddo
	x(jiR)=t+b(jiR)
    call MPI_ALLREDUCE(x(ji:jf),b(ji:jf),Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	x(ji:jf)=b(ji:jf)
	do j=ji,jf
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=ji,j-1
	    t=t-conjg(A(k,j))*x(k)
      enddo
	  x(j)=t
	enddo
  enddo
  mod_n_Glob_NumOfProcs=mod(n,Glob_NumOfProcs)
  if (mod_n_Glob_NumOfProcs>0) then
    jim=RowsPerProc*Glob_NumOfProcs
	ji=jim+1
	jiR=ji+Glob_ProcID
	if (jiR<=n) then
	  t=ZERO
	  !DEC$ UNROLL (4)
      do k=1,jim
        t=t-conjg(A(k,jiR))*x(k)
	  enddo
	  x(jiR)=t+b(jiR)
    endif
    call MPI_ALLREDUCE(x(ji:n),b(ji:n),mod_n_Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    x(ji:n)=b(ji:n)
	do j=ji,n
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=ji,j-1
	    t=t-conjg(A(k,j))*x(k)
      enddo
	  x(j)=t
	enddo
  endif
  !Solution of D*LH*x=y  
  b(1:n)=x(1:n)*invD(1:n)   
  x(1:n)=ZERO
  do i=1,RowsPerProc
    jf=n-i*Glob_NumOfProcs+1
	ji=jf-1+Glob_NumOfProcs
	jiR=ji+1+Glob_ProcID
	if (jiR<=n) then
 	  t=x(jiR)
 	  !DEC$ UNROLL (4)
      do k=ji,1,-1
         x(k)=x(k)-t*A(k,jiR)
	  enddo
      call MPI_ALLREDUCE(x(jf:ji),b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs),Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      x(jf:ji)=b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
	endif
    x(jf:ji)=x(jf:ji)+b(jf:ji)
	do j=ji,jf,-1
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=j-1,jf,-1
	     x(k)=x(k)-t*A(k,j)
	  enddo
	enddo
  enddo
  if (mod_n_Glob_NumOfProcs>0) then
	ji=mod_n_Glob_NumOfProcs
	jiR=ji+1+Glob_ProcID
	if (jiR<=n) then
 	  t=x(jiR)
 	  !DEC$ UNROLL (4)
      do k=ji,1,-1
         x(k)=x(k)-t*A(k,jiR)
	  enddo
      call MPI_ALLREDUCE(x(1:ji),b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs),mod_n_Glob_NumOfProcs,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      x(1:ji)=b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
	endif
	x(1:ji)=x(1:ji)+b(1:ji)
	do j=ji,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
	  do k=j-1,1,-1
	     x(k)=x(k)-t*A(k,j)
	  enddo
	enddo
  endif
endif

end subroutine LDLHS



subroutine MTMVL(n,A,nA,x,y,w)
!Subroutine MTMVL computes the product of real symmetric matrix
!A and vector x using only the lower triangle of A. The upper
!triangle is not referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A
!           (its lower triangle, including the diagonal);
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of lenght n)
!  Output parameters :
!       y - An array (vector) containing the result; 

integer      n,nA,j,jm,ji,jf,p,i
real(dprec)  A(nA,n),x(n),y(n),w(n),t

 if (Glob_MTMVL_PMode==0) then
  if (Glob_MTMVL_UseBLAS==0) then  
    !We use the property A*x=(R+L)*x = (x^T*R^T)^T + (x^T*L^T)^T
    !where L is the lower triangle of A (without the diagonal) and
    !R is the upper triangle (including the diagonal).
    !Computing (x^T*R^T)^T
    do j=1,n
      t=ZERO
      !DEC$ UNROLL (4)
      do i=j,n
        t=t+A(i,j)*x(i)
      enddo
      y(j)=t
    enddo  
    !Computing (x^T*L^T)^T
    do j=n,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
      do i=j+1,n
        y(i)=y(i)+t*A(i,j)
      enddo
    enddo
  else
    !call BLAS routine DSYMV
    call DSYMV('L',n,ONE,A,nA,x,1,ZERO,y,1)
  endif
else
  w(1:n)=ZERO
  p=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) p=p+1
  !Computing (x^T*R^T)^T
  ji=1+p*Glob_ProcID
  jf=min(p*(Glob_ProcID+1),n)
  do j=ji,jf
    t=ZERO
    !DEC$ UNROLL (4)
    do i=j,n
      t=t+A(i,j)*x(i)
    enddo
    w(j)=t
  enddo
  !Computing (x^T*L^T)^T
  ji=1+p*(Glob_NumOfProcs-Glob_ProcID-1)
  jf=min(p*(Glob_NumOfProcs-Glob_ProcID),n)
  do j=jf,ji,-1
    t=x(j)
    !DEC$ UNROLL (4)
    do i=j+1,n
      w(i)=w(i)+t*A(i,j)
    enddo
  enddo
  call MPI_ALLREDUCE(w,y,n,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

end subroutine MTMVL


subroutine MHMVL(n,A,nA,x,y,w)
!Subroutine MHMVL computes the product of complex hermitian matrix
!A and vector x using only the lower triangle of A. The upper
!triangle is not referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A
!           (its lower triangle, including the diagonal);
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of lenght n)
!  Output parameters :
!       y - An array (vector) containing the result; 

integer         n,nA,j,jm,ji,jf,p,i
complex(dprec)  A(nA,n),x(n),y(n),w(n),t

if (Glob_MHMVL_PMode==0) then
  if (Glob_MHMVL_UseBLAS==0) then
    !We use the property A*x=(R+L)*x = (x^H*R^H)^H + (x^H*L^H)^H
    !where L is the lower triangle of A (without the diagonal) and
    !R is the upper triangle (including the diagonal).   
    !Computing (x^H*R^H)^H
    do j=1,n
      t=cmplx(ZERO,ZERO,dprec)
      !DEC$ UNROLL (4)
      do i=j,n
        t=t+conjg(A(i,j))*x(i)
      enddo
      y(j)=t    
    enddo
    !Computing (x^H*L^H)^H
    do j=n,1,-1
      t=x(j)
      !DEC$ UNROLL (4)
      do i=j+1,n
        y(i)=y(i)+t*A(i,j)
      enddo
    enddo
  else
    !call BLAS routine ZHEMV
    call ZHEMV('L',n,cmplx(ONE,ZERO,dprec),A,nA,x,1,cmplx(ZERO,ZERO,dprec),y,1)  
  endif
else
  w(1:n)=ZERO
  p=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) p=p+1
  !Computing (x^H*R^H)^H
  ji=1+p*Glob_ProcID
  jf=min(p*(Glob_ProcID+1),n)
  do j=ji,jf
    t=cmplx(ZERO,ZERO,dprec)
    !DEC$ UNROLL (4)
    do i=j,n
      t=t+conjg(A(i,j))*x(i)
    enddo
    w(j)=t
  enddo
  !Computing (x^H*L^H)^H
  ji=1+p*(Glob_NumOfProcs-Glob_ProcID-1)
  jf=min(p*(Glob_NumOfProcs-Glob_ProcID),n)
  do j=jf,ji,-1
    t=x(j)
    !DEC$ UNROLL (4)
    do i=j+1,n
      w(i)=w(i)+t*A(i,j)
    enddo
  enddo
  call MPI_ALLREDUCE(w,y,n,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

end subroutine MHMVL



subroutine MTMV(n,A,nA,x,y,w)
!Subroutine MTMV computes the product of transposed 
!real matrix A and vector x. 
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A;
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of lenght n)
!  Output parameters :
!       y - An array (vector) containing the result, y=A^T*x  

integer      n,nA,ib,ie,jb,je,j,i,k
real(dprec)  A(nA,n),x(n),y(n),w(n),t
integer      n2,p,q

if (Glob_MTMV_PMode==0) then
  if (Glob_MTMV_UseBLAS==0) then
    do i=1,n
      t=ZERO
      !DEC$ UNROLL (4)
      do j=1,n
        t=t+A(j,i)*x(j)
      enddo
      y(i)=t
    enddo
  else
    !call BLAS routine DGEMV
    call DGEMV('T',n,n,ONE,A,nA,x,1,ZERO,y,1)
  endif  
else
  !Each process is assigned to deal with some part of 
  !matrix A. This part includes all matrix elements between
  !ib,jb and  ie,je (i and j stand for row and index 
  !respectively). The word "between" means in natural order, 
  !A11, A12, ..., A1n, A21, A22, ... Since we actually need the 
  !transposed matrix, the access to matrix elements will therefore 
  !be by rows, which suitable for Fortran type of storage. The 
  !actual number of elements assigned to a particular process 
  !is k. The algorithm divides the labor between processes as 
  !uniformly as possible.
  w(1:n)=ZERO
  n2=n*n
  p=n2/Glob_NumOfProcs
  if (mod(n2,Glob_NumOfProcs)/=0) p=p+1
  q=Glob_ProcID*p
  if (q+p>n2) then 
    k=max(n2-q,0)
  else
    k=p
  endif
  ib=q/n+1
  ie=min((q+p-1)/n+1,n)
  jb=mod(q,n)+1
  je=mod(q+k-1,n)+1
  if (ie>ib) then
    t=ZERO
    !DEC$ UNROLL (4)
    do j=jb,n
      t=t+A(j,ib)*x(j)
    enddo
    w(ib)=t
    do i=ib+1,ie-1
      t=ZERO
      !DEC$ UNROLL (4)
      do j=1,n
        t=t+A(j,i)*x(j)
      enddo
      w(i)=t
    enddo
    t=ZERO
    !DEC$ UNROLL (4)
    do j=1,je
      t=t+A(j,i)*x(j)
    enddo
    w(ie)=t
  else
    if (ie==ib) then
      t=ZERO
      !DEC$ UNROLL (4)
      do j=jb,je
        t=t+A(j,ib)*x(j)
      enddo
      w(ib)=t
    endif  
  endif
  call MPI_ALLREDUCE(w,y,n,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end subroutine MTMV


subroutine MHMV(n,A,nA,x,y,w)
!Subroutine MHMV computes the product of conjugated and transposed 
!complex matrix A and vector x. 
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A;
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of lenght n)
!  Output parameters :
!       y - An array (vector) containing the result, y=A^H*x  

integer         n,nA,ib,ie,jb,je,j,i,k
complex(dprec)  A(nA,n),x(n),y(n),w(n),t
integer         n2,p,q

if (Glob_MHMV_PMode==0) then
  if (Glob_MHMV_UseBLAS==0) then
    do i=1,n
      t=cmplx(ZERO,ZERO,dprec)
      !DEC$ UNROLL (4)
      do j=1,n
        t=t+conjg(A(j,i))*x(j)
      enddo
      y(i)=t
    enddo
  else
    !call BLAS routine ZGEMV
    call ZGEMV('C',n,n,cmplx(ONE,ZERO,dprec),A,nA,x,1,cmplx(ZERO,ZERO,dprec),y,1) 
  endif
else
  !Each process is assigned to deal with some part of 
  !matrix A. This part includes all matrix elements between
  !ib,jb and  ie,je (i and j stand for row and index 
  !respectively). The word "between" means in natural order, 
  !A11, A12, ..., A1n, A21, A22, ...  Since we actually need 
  !the transposed (and conjugated) matrix, the access to 
  !matrix elements will therefore be by rows, which suitable
  !for Fortran type of storage. The actual number of elements 
  !assigned to a particular process is k. The algorithm
  !divides the labor between processes as uniformly as possible.
  w(1:n)=cmplx(ZERO,ZERO,dprec)
  n2=n*n
  p=n2/Glob_NumOfProcs
  if (mod(n2,Glob_NumOfProcs)/=0) p=p+1
  q=Glob_ProcID*p
  if (q+p>n2) then 
    k=max(n2-q,0)
  else
    k=p
  endif
  ib=q/n+1
  ie=min((q+p-1)/n+1,n)
  jb=mod(q,n)+1
  je=mod(q+k-1,n)+1
  if (ie>ib) then
    t=cmplx(ZERO,ZERO,dprec)
    !DEC$ UNROLL (4)
    do j=jb,n
      t=t+conjg(A(j,ib))*x(j)
    enddo
    w(ib)=t
    do i=ib+1,ie-1
      t=cmplx(ZERO,ZERO,dprec)
      !DEC$ UNROLL (4)
      do j=1,n
        t=t+conjg(A(j,i))*x(j)
      enddo
      w(i)=t
    enddo
    t=cmplx(ZERO,ZERO,dprec)
    !DEC$ UNROLL (4)
    do j=1,je
      t=t+conjg(A(j,i))*x(j)
    enddo
    w(ie)=t
  else
    if (ie==ib) then
      t=cmplx(ZERO,ZERO,dprec)
      !DEC$ UNROLL (4)
      do j=jb,je
        t=t+conjg(A(j,ib))*x(j)
      enddo
      w(ib)=t
    endif  
  endif
  call MPI_ALLREDUCE(w,y,n,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)  
endif
end subroutine MHMV


function VMMTMV(n,A,nA,x)
!Function VMMTMV computes the product x^T*A*x, where x is a real vector,
!x^T is transposed x, and A is a real symmetric matrix.
!Only the lower triangle of A (including the diagonal) is referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A (only the lower 
!           triangle is needed);
!      nA - The leading dimension of A;
!       x - An array containing vector x;

!Arguments
integer        n,nA
real(dprec)    A(nA,n),x(n)
real(dprec)    VMMTMV
!Local variables  
integer      j,k,q,p
real(dprec)  s,sum

if (Glob_VMMTMV_PMode==0) then
  VMMTMV=ZERO
  do k=1,n
    VMMTMV=VMMTMV+x(k)*x(k)*A(k,k)
	s=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  s=s+x(j)*A(j,k)
 	enddo
    VMMTMV=VMMTMV+TWO*s*x(k)
  enddo
else
  sum=ZERO
  q=n/(2*Glob_NumOfProcs)
  do p=1,2*q,2
    k=(p-1)*Glob_NumOfProcs+Glob_ProcID+1
    sum=sum+x(k)*x(k)*A(k,k)
	s=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  s=s+x(j)*A(j,k)
	enddo
	sum=sum+TWO*s*x(k)
	k=(p+1)*Glob_NumOfProcs-Glob_ProcID   
    sum=sum+x(k)*x(k)*A(k,k)
	s=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  s=s+x(j)*A(j,k)
	enddo
	sum=sum+TWO*s*x(k)
  enddo  
  p=Glob_ProcID
  do k=2*q*Glob_NumOfProcs+1,n
	if (p==0) then
      sum=sum+x(k)*x(k)*A(k,k)
	  p=p+Glob_NumOfProcs
	endif
	s=ZERO
    do j=k+p,n,Glob_NumOfProcs
	  s=s+x(j)*A(j,k)	
	enddo
	sum=sum+TWO*s*x(k)
	p=j-n-1 
  enddo
  call MPI_ALLREDUCE(sum,VMMTMV,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

end function VMMTMV


function VMMHMV(n,A,nA,x)
!Function VMMHMV computes the product x^H*A*x, where x is a complex vector,
!x^H is conjugated and transposed x, and A is a complex hermitian matrix.
!Only the lower triangle of A (including the diagonal) is referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A (only the lower 
!           triangle is needed);
!      nA - The leading dimension of A;
!       x - An array containing vector x;

!Arguments
integer        n,nA
complex(dprec) A(nA,n),x(n)
real(dprec)    VMMHMV
!Local variables  
integer      j,k,q,p
real(dprec)  sr,si,sum

if (Glob_VMMHMV_PMode==0) then
  VMMHMV=ZERO
  do k=1,n
    VMMHMV=VMMHMV+(real(x(k),dprec)*real(x(k),dprec)+imag(x(k))*imag(x(k)))*real(A(k,k),dprec)
	sr=ZERO
	si=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  sr=sr+(real(x(j),dprec)*real(A(j,k),dprec)+imag(x(j))*imag(A(j,k)))
	  si=si+(imag(x(j))*real(A(j,k),dprec)-real(x(j),dprec)*imag(A(j,k)))
 	enddo
    VMMHMV=VMMHMV+TWO*(real(x(k),dprec)*sr+imag(x(k))*si)
  enddo
else
  sum=ZERO
  q=n/(2*Glob_NumOfProcs)
  do p=1,2*q,2
    k=(p-1)*Glob_NumOfProcs+Glob_ProcID+1
    sum=sum+(real(x(k),dprec)*real(x(k),dprec)+imag(x(k))*imag(x(k)))*real(A(k,k),dprec)
	sr=ZERO
	si=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  sr=sr+(real(x(j),dprec)*real(A(j,k),dprec)+imag(x(j))*imag(A(j,k)))
	  si=si+(imag(x(j))*real(A(j,k),dprec)-real(x(j),dprec)*imag(A(j,k)))
	enddo
	sum=sum+TWO*(real(x(k),dprec)*sr+imag(x(k))*si)
	k=(p+1)*Glob_NumOfProcs-Glob_ProcID   
    sum=sum+(real(x(k),dprec)*real(x(k),dprec)+imag(x(k))*imag(x(k)))*real(A(k,k),dprec)
	sr=ZERO
	si=ZERO
	!DEC$ UNROLL (4)
    do j=k+1,n
	  sr=sr+(real(x(j),dprec)*real(A(j,k),dprec)+imag(x(j))*imag(A(j,k)))
	  si=si+(imag(x(j))*real(A(j,k),dprec)-real(x(j),dprec)*imag(A(j,k)))
	enddo
	sum=sum+TWO*(real(x(k),dprec)*sr+imag(x(k))*si)
  enddo  
  p=Glob_ProcID
  do k=2*q*Glob_NumOfProcs+1,n
	if (p==0) then
      sum=sum+(real(x(k),dprec)*real(x(k),dprec)+imag(x(k))*imag(x(k)))*real(A(k,k),dprec)
	  p=p+Glob_NumOfProcs
	endif
	sr=ZERO
	si=ZERO
    do j=k+p,n,Glob_NumOfProcs
	  sr=sr+(real(x(j),dprec)*real(A(j,k),dprec)+imag(x(j))*imag(A(j,k)))
	  si=si+(imag(x(j))*real(A(j,k),dprec)-real(x(j),dprec)*imag(A(j,k)))	
	enddo
	sum=sum+TWO*(real(x(k),dprec)*sr+imag(x(k))*si)
	p=j-n-1 
  enddo
  call MPI_ALLREDUCE(sum,VMMHMV,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif

end function VMMHMV


function RMaxAbsEl(n,x)
!Function RMaxAbsEl returns the magnitude of the largest by magnitude element 
!among the first n elements of real array x
integer        n,j,ji,jf,k
real(dprec)    x(n)
real(dprec)    RMaxAbsEl,MaxAE
MaxAE=ZERO
if (Glob_RMaxAbsEl_PMode==0) then
  do j=1,n
    if (abs(x(j))>MaxAE) MaxAE=abs(x(j))
  enddo
  RMaxAbsEl=MaxAE
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  do j=ji,jf
    if (abs(x(j))>MaxAE) MaxAE=abs(x(j))
  enddo
  call MPI_ALLREDUCE(MaxAE,RMaxAbsEl,1,MPI_DPREC,MPI_MAX,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function RMaxAbsEl


function CMaxAbsReOrIm(n,x)
!Function CMaxAbsEl returns the magnitude of the largest real or imaginary part 
!among the first n elements of complex array x
integer  n,j,ji,jf,k,maxj
complex(dprec)  x(n)
real(dprec)     CMaxAbsReOrIm,MaxAE,t
MaxAE=ZERO
if (Glob_CMaxAbsReOrIm_PMode==0) then
  do j=1,n
    t=abs(real(x(j),dprec))
    if (t>MaxAE) MaxAE=t
    t=abs(imag(x(j)))
    if (t>MaxAE) MaxAE=t
  enddo
  CMaxAbsReOrIm=MaxAE
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  do j=ji,jf
    t=abs(real(x(j),dprec))
    if (t>MaxAE) MaxAE=t
    t=abs(imag(x(j)))
    if (t>MaxAE) MaxAE=t
  enddo
  call MPI_ALLREDUCE(MaxAE,CMaxAbsReOrIm,1,MPI_DPREC,MPI_MAX,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function CMaxAbsReOrIm


function RDotProd(n,x,y)
!Function RDotProd computes the dot product x^{T}y,
!for n-component real vectors x and y.
integer           n,k,ji,jf,j
real(dprec)    RDotProd,a
real(dprec)    x(n),y(n)
real(dprec)    DDOT
if (Glob_RDotProd_PMode==0) then
  if (Glob_RDotProd_UseBLAS==0) then
    RDotProd=ZERO
    !DEC$ UNROLL (4)
    do j=1,n
      RDotProd=RDotProd+x(j)*y(j)
    enddo
  else
    !call BLAS function DDOT
     RDotProd=DDOT(n,x,1,y,1)
  endif  
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  a=ZERO
  !DEC$ UNROLL (4)
  do j=ji,jf
    a=a+x(j)*y(j)
  enddo
  call MPI_ALLREDUCE(a,RDotProd,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function RDotProd


function CDotProd(n,x,y)
!Function CDotProd computes the dot product x^{H}y,
!for n-component complex vectors x and y.
integer           n,k,ji,jf,j
complex(dprec)    CDotProd,a
complex(dprec)    x(n),y(n)
complex(dprec)    ZDOTC
if (Glob_CDotProd_PMode==0) then
  if (Glob_CDotProd_UseBLAS==0) then
    CDotProd=cmplx(ZERO,ZERO,dprec)
    !DEC$ UNROLL (4)
    do j=1,n
      CDotProd=CDotProd+cmplx(real(x(j),dprec)*real(y(j),dprec)+imag(x(j))*imag(y(j)),real(x(j),dprec)*imag(y(j))-imag(x(j))*real(y(j),dprec),dprec)
    enddo
  else
    !call BLAS function ZDOTC
    CDotProd=ZDOTC(n,x,1,y,1)
  endif
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  a=cmplx(ZERO,ZERO,dprec)
  !DEC$ UNROLL (4)
  do j=ji,jf
    a=a+cmplx(real(x(j),dprec)*real(y(j),dprec)+imag(x(j))*imag(y(j)),real(x(j),dprec)*imag(y(j))-imag(x(j))*real(y(j),dprec),dprec)
  enddo
  call MPI_ALLREDUCE(a,CDotProd,1,MPI_DCOMP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function CDotProd


function RDotProdItself(n,x)
!Function RDotProdItself computes dot product x^{T}x,
!for an n-component real vector x.
integer           n,k,ji,jf,j
real(dprec)       RDotProdItself,t
real(dprec)       x(n)
real(dprec)       DDOT
if (Glob_RDotProdItself_PMode==0) then
  if (Glob_RDotProdItself_UseBLAS==0) then
    RDotProdItself=ZERO
    !DEC$ UNROLL (4)
    do j=1,n
      RDotProdItself=RDotProdItself+x(j)*x(j)
    enddo
  else
    !call BLAS function DDOT
    RDotProdItself=DDOT(n,x,1,x,1)
  endif  
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  t=ZERO
  !DEC$ UNROLL (4)
  do j=ji,jf
    t=t+x(j)*x(j) 
  enddo
  call MPI_ALLREDUCE(t,RDotProdItself,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function RDotProdItself


function CDotProdItself(n,x)
!Function CDotProdItself computes dot product x^{H}x,
!for an n-component complex vector x.
integer           n,k,ji,jf,j
real(dprec)       CDotProdItself,t
complex(dprec)    x(n)
complex(dprec)    ZDOTC
if (Glob_CDotProdItself_PMode==0) then
  if (Glob_CDotProdItself_UseBLAS==0) then  
    CDotProdItself=ZERO
    !DEC$ UNROLL (4)
    do j=1,n
      CDotProdItself=CDotProdItself+real(x(j),dprec)*real(x(j),dprec)+imag(x(j))*imag(x(j))
    enddo
  else
    !call BLAS function ZDOTC
    CDotProdItself=ZDOTC(n,x,1,x,1)
  endif  
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  t=ZERO
  !DEC$ UNROLL (4)
  do j=ji,jf
    t=t+real(x(j),dprec)*real(x(j),dprec)+imag(x(j))*imag(x(j))    
  enddo
  call MPI_ALLREDUCE(t,CDotProdItself,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
endif
end function CDotProdItself



function RDotProdQuotient(n,x,y)
!Function RDotProdQuotient computes the quotion (x^{T}y)/(y^{T}y),
!where x and y are n-component real vectors.
integer        n,k,ji,jf,j
real(dprec)    RDotProdQuotient,a,b,t(2),tt(2)
real(dprec)    x(n),y(n)
if (Glob_RDotProdQuotient_PMode==0) then
  a=ZERO
  b=ZERO
  !DEC$ UNROLL (4)
  do j=1,n
    a=a+x(j)*y(j)
    b=b+y(j)*y(j)
  enddo 
  RDotProdQuotient=a/b
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  t(1)=ZERO
  t(2)=ZERO
  !DEC$ UNROLL (4)
  do j=ji,jf
    t(1)=t(1)+x(j)*y(j)
    t(2)=t(2)+y(j)*y(j)
  enddo
  call MPI_ALLREDUCE(t,tt,2,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  RDotProdQuotient=tt(1)/tt(2)
endif
end function RDotProdQuotient



function CDotProdQuotient(n,x,y)
!Function CDotProdQuotient computes the quotion (x^{H}y)/(y^{H}y),
!where x and y are n-component complex vectors.
integer        n,k,ji,jf,j
complex(dprec) CDotProdQuotient
complex(dprec) x(n),y(n)
complex(dprec) a
real(dprec)    r,t(3),tt(3)
if (Glob_CDotProdQuotient_PMode==0) then
  a=cmplx(ZERO,ZERO,dprec)
  r=ZERO
  !DEC$ UNROLL (4)
  do j=1,n
    a=a+conjg(x(j))*y(j)
    r=r+real(y(j),dprec)*real(y(j),dprec)+imag(y(j))*imag(y(j))
  enddo 
  CDotProdQuotient=cmplx(real(a,dprec)/r,imag(a)/r,dprec)
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  a=cmplx(ZERO,ZERO,dprec)
  t(1)=ZERO
  !DEC$ UNROLL (4)
  do j=ji,jf
    a=a+conjg(x(j))*y(j)
    t(1)=t(1)+real(y(j),dprec)*real(y(j),dprec)+imag(y(j))*imag(y(j))
  enddo
  t(2)=real(a,dprec)
  t(3)=imag(a)
  call MPI_ALLREDUCE(t,tt,3,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  CDotProdQuotient=cmplx(tt(2)/tt(1),tt(3)/tt(1),dprec)
endif
end function CDotProdQuotient



subroutine RVScale(n,x,alpha)
!Subroutine RVScale scales real vector x by real constant alpha
integer      n,j
real(dprec)  x(n)
real(dprec)  alpha

if (Glob_RVScale_UseBLAS==0) then
  !DEC$ UNROLL (4)
  do j=1,n
    x(j)=alpha*x(j)
  enddo  
else
  !call BLAS routine DSCAL
  call DSCAL(n,alpha,x,1)
endif

end subroutine RVScale


subroutine CVScale(n,x,alpha)
!Subroutine RVScale scales complex vector x by complex constant alpha
integer         n,j
complex(dprec)  x(n)
complex(dprec)  alpha
real(dprec)     ralpha

if (Glob_CVScale_UseBLAS==0) then
  if (imag(alpha)==ZERO) then
    ralpha=real(alpha,dprec)
    !DEC$ UNROLL (4)
    do j=1,n
      x(j)=cmplx(ralpha*real(x(j),dprec),ralpha*imag(x(j)),dprec)
    enddo
  else
    !DEC$ UNROLL (4)
    do j=1,n
      x(j)=x(j)*alpha
    enddo    
  endif  
else
  !call BLAS routine ZSCAL
  call ZSCAL(n,alpha,x,1)
endif

end subroutine CVScale



function RVDiffEucNorm(n,x,alpha,y)
!Function RVDiffEucNorm returns the Euclidian norm of the difference x-alpha*y, 
!where y and x are two real vectors and alpha is a real scalar.
real(dprec)  RVDiffEucNorm
integer      n,j,k,ji,jf,q
real(dprec)  x(n),y(n),alpha,t

if (Glob_RVDiffEucNorm_PMode==0) then
  RVDiffEucNorm=ZERO    
  !DEC$ UNROLL (4)
  do j=1,n    
      RVDiffEucNorm=RVDiffEucNorm+(x(j)-alpha*y(j))*(x(j)-alpha*y(j))
  enddo       
  RVDiffEucNorm=sqrt(RVDiffEucNorm)
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  t=ZERO 
  !DEC$ UNROLL (4)
  do j=ji,jf
      t=t+(x(j)-alpha*y(j))*(x(j)-alpha*y(j)) 
  enddo
  call MPI_ALLREDUCE(t,RVDiffEucNorm,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  RVDiffEucNorm=sqrt(RVDiffEucNorm)
endif

end function RVDiffEucNorm


function CVDiffEucNorm(n,x,alpha,y)
!Function CVDiffEucNorm returns the Euclidian norm of the difference x-alpha*y, 
!where y and x are two complex vectors and alpha is a complex scalar.
real(dprec)     CVDiffEucNorm
integer         n,j,k,ji,jf,q
complex(dprec)  x(n),y(n),alpha,s
real(dprec)     t

if (Glob_CVDiffEucNorm_PMode==0) then
  CVDiffEucNorm=ZERO    
  !DEC$ UNROLL (4)
  do j=1,n    
    s=x(j)-alpha*y(j)
    CVDiffEucNorm=CVDiffEucNorm+real(s,dprec)*real(s,dprec)+imag(s)*imag(s)
  enddo       
  CVDiffEucNorm=sqrt(CVDiffEucNorm)
else
  k=n/Glob_NumOfProcs
  if (mod(n,Glob_NumOfProcs)/=0) k=k+1
  ji=1+k*Glob_ProcID
  jf=min(k*(Glob_ProcID+1),n)
  t=cmplx(ZERO,ZERO,dprec)
  !DEC$ UNROLL (4)
  do j=ji,jf
    s=x(j)-alpha*y(j)
    t=t+real(s,dprec)*real(s,dprec)+imag(s)*imag(s)
  enddo
  call MPI_ALLREDUCE(t,CVDiffEucNorm,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  CVDiffEucNorm=sqrt(CVDiffEucNorm)
endif

end function CVDiffEucNorm



subroutine GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
    lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
!Subroutine GSEPIIS finds a nondegenerate eigenvalue and the
!corresponding eigenvector of the generalized real symmetric eigenvalue
!problem A*x=lambda*B*x with symmetric matrices A and B using the
!inverse iteration method. The subroutine makes use of a 
!known factorization of A-apprlambda*B submatrix of size k-1
!in the form M=A-apprlambda*B=L*D*LT which is done by 
!routine LDLTF. The desired eigenvalue must not have any other
!pathologically close eigenvalues (i.e. with relative separation 
!of the order MachineEpsilon and less), otherwise the solution procees 
!cannot separate them from each other and they will be interpreted 
!as one degenerate eigenvalue. In the last case the eigenvector will 
!lie in the subspace spanned by the corresponding eigenvectors 
!and the subroutine can generate ErrorCode=2 upon exit. In 
!the case of degenerate eigenvalues the accuracy is not guaranteed,
!however if a group of eigenvalues is much closer 
!to apprlambda than other eigenvalues then after sufficient 
!number of iterations desirable eigenvector will be found, 
!even though there may be generated ErrorCode=2.
!  Input parameters :
!     k-1 - The size of M=A-apprlambda*B submatrix, whose
!           factorization is known (k>=1);
!       n - The size of matrices A and B;
!       M - A two-dimensional array containing known submatrix LH 
!           (of size k-1) in the upper triangle (excluding the 
!           diagonal, the latter is assumed to consist of ones)
!           and submatrix A-apprlambda*B (of the size k-1) in 
!           the lower triangle (including the diagonal); 
!      nM - The leading dimension of M;
!    invD - An array containing the inversed values of the 
!           diagonal matrix D (of size k-1);
!       B - A two-dimensional array containing matrix B;
!      nB - The leading dimension of B;
!apprlambda - An approximation to the desired eigenvalue. To
!           get a converging process, it is necessary to have 
!           the distance between apprlambda and the desired 
!           eigenvalue smaller than the distance between apprlambda
!           and any other eigenvalue. The closer apprlabda to
!           the desired eigenvalue the faster the iterative 
!           process will converge. However, apprlambda cannot be
!           exactly equal to the real eigenvalue as in this case the 
!           matrices the routine deals with are going to be be singular;
!       v - An array (at least of size n) where an approximation
!           to the desired eigenvector should be placed. Having
!           a good approximation helps to speed up the iterative process.
!           If such an approximation is not known, it is quite 
!           acceptable to use a random vector (but not a zero vector!). 
!           This array is also used as a work array. Thus, it is 
!           destroyed upon exit;
!       w - A work array (at least of size n);
!     Tol - The desired relative accuracy of the calculations (i.e. the
!           relative accuracy of the eigevector and the eigenvalue).
!           If this value is set to be small and negative then the 
!           routine finds as accurate solution as possible, but not less
!           accurate than abs(Tol). In the case of negative Tol
!           the calculations require at least one more iteration (but potentially
!           can result in many more iterations, up to the iteration limit, MaxIter)
!           and the total computation time may increase substantially.
! MaxIter - The maximum number of iterations allowed;
!SpecifNorm - The parameter that determines how the eivenvector is 
!           normalized upon exit. If SpecifNorm=0 then it is normalized 
!           in such a way that xT*B*x=1. If SpecifNorm=1 then xT*x=1. 
!           Otherwise the eigenvector is normalized in such a way that
!           the largest by magnitude component is one; 
!  Output parameters:
!       M - A two-dinensional array whose upper triangle (excluding
!           the diagonal) contains the updated matrix LT (of
!           size n) and the lower triangle (including the 
!           diagonal) contains the lower triangle of matrix
!           A-apprlambda*B (of size n); 
!    invD - An array containing the inversed values of
!           the diagonal matrix D (of size n);
!  lambda - The eigenvalue found;
!       x - The eigenvector found;
!  RelAcc - A rough estimate of the relative accuracy reached.
!           It may be somewhat inaccurate. A more accurate estimate
!           would require additional calculations and for the
!           sake of efficiency was not implemented.
! NumIter - The number of iterations performed;
!ErrorCode - An error flag. If ErrorCode=0 then the routine
!           finished succesfully. There may be the following
!           errors on exit:
!           ErrorCode=1 - matrix A-apprlambda*B is singular 
!           or almost singular;
!           ErrorCode=2 - the process did not converge to the 
!           accuracy Tol after the maximum number of 
!           iterations was performed;           

!Arguments
integer         k,n,nM,nB,MaxIter,SpecifNorm,NumIter,ErrorCode
real(dprec)     M(nM,n),invD(n),B(nB,n),v(n),w(n),x(n)
real(dprec)     apprlambda,Tol,lambda,RelAcc
!Local variables
real(dprec)     NormOfDiff,NormOfDiffPrev,t1,t2,sqrtn
real(dprec)     tc
logical         notconverged

ErrorCode=0
RelAcc=huge(RelAcc)
NormOfDiffPrev=huge(NormOfDiffPrev)
NumIter=0
!Updating M=L*D*LH factorization up to size n using routine LDLHF
call LDLTF(k,n,M,nM,invD,w,ErrorCode)
if (ErrorCode>0) return
sqrtn=sqrt(n*ONE)
!Do inverse iterations until the process converges
!with relative accuracy Tol, or until the number
!of iterations exceeds the limit
notconverged=.true.
do while ((notconverged).AND.(NumIter<MaxIter))
  if (NumIter/=0) v(1:n)=x(1:n)
  call MTMV(n,B,nB,v,w,x)
  call LDLTS(n,M,nM,invD,w,x)
  t1=RMaxAbsEl(n,x)
  call RVScale(n,x,1/t1)  
  !tc=(x^{T}v)/(v^{T}v)
  tc=RDotProdQuotient(n,x,v) 
  !NormOfDiff=sqrt(n)||x-tc*v||
  NormOfDiff=sqrtn*RVDiffEucNorm(n,x,tc,v)
  if (Tol>ZERO) then
    if (NormOfDiff<=Tol) notconverged=.false.
  else
    if ((NormOfDiff>NormOfDiffPrev).AND.(NormOfDiff<=abs(Tol))) notconverged=.false.
    NormOfDiffPrev=NormOfDiff
  endif  
  NumIter=NumIter+1
enddo
RelAcc=NormOfDiff
if (notconverged) ErrorCode=2
!Compute x^{T}Bx
t1=VMMTMV(n,B,nB,x)
!Compute x^{T}Mx
t2=VMMTMV(n,M,nM,x)
!Rayleigh quotient
lambda=(t2/t1)+apprlambda  
select case	(SpecifNorm)
  case (0) !x^{H}Bx=1
    call RVScale(n,x,1/sqrt(t1))
  case (1) !x^{H}x=1 
    t1=RDotProdItself(n,x)
    call RVScale(n,x,1/sqrt(t1))
endselect

end subroutine GSEPIIS



subroutine GHEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
    lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
!Subroutine GHEPIIS finds a nondegenerate eigenvalue and the
!corresponding eigenvector of the generalized complex hermitian eigenvalue
!problem A*x=lambda*B*x with hermitian matrices A and B using the
!inverse iteration method. The subroutine makes use of a 
!known factorization of A-apprlambda*B submatrix of size k-1
!in the form A-apprlambda*B=L*D*LH which is done by 
!routine LDLHF. The desired eigenvalue must not have any other
!pathologically close eigenvalues (i.e. with relative separation 
!of the order MachineEpsilon and less), otherwise the solution procees 
!cannot separate them from each other and they will be interpreted 
!as one degenerate eigenvalue. In the last case the eigenvector will 
!lie in the subspace spanned by the corresponding eigenvectors 
!and the subroutine can generate ErrorCode=2 upon exit. In 
!the case of degenerate eigenvalues the accuracy is not guaranteed,
!however if a group of eigenvalues is much closer 
!to apprlambda than other eigenvalues then after sufficient 
!number of iterations desirable eigenvector will be found, 
!even though there may be generated ErrorCode=2.
!  Input parameters :
!     k-1 - The size of A-apprlambda*B submatrix, whose
!           factorization is known (k>=1);
!       n - The size of matrices A and B;
!       M - A two-dimensional array containing known submatrix LH 
!           (of size k-1) in the upper triangle (excluding the 
!           diagonal, the latter is assumed to consist of ones)
!           and submatrix A-apprlambda*B (of the size k-1) in 
!           the lower triangle (including the diagonal); 
!      nM - The leading dimension of M;
!    invD - An array containing the inversed values of the 
!           diagonal matrix D (of size k-1);
!       B - A two-dimensional array containing matrix B;
!      nB - The leading dimension of B;
!apprlambda - An approximation to the desired eigenvalue. To
!           get a converging process, it is necessary to have 
!           the distance between apprlambda and the desired 
!           eigenvalue smaller than the distance between apprlambda
!           and any other eigenvalue. The closer apprlabda to
!           the desired eigenvalue the faster the iterative 
!           process will converge. However, apprlambda cannot be
!           exactly equal to the real eigenvalue as in this case the 
!           matrices the routine deals with are going to be be singular;
!       v - An array (at least of size n) where an approximation
!           to the desired eigenvector should be placed. Having
!           a good approximation helps to speed up the iterative process.
!           If such an approximation is not known, it is quite 
!           acceptable to use a random vector (but not a zero vector!). 
!           This array is also used as a work array. Thus, it is 
!           destroyed upon exit;
!       w - A work array (at least of size n);
!     Tol - The desired relative accuracy of the calculations (i.e. the
!           relative accuracy of the eigevector and the eigenvalue).
!           If this value is set to be small and negative then the 
!           routine finds as accurate solution as possible, but not less
!           accurate than abs(Tol). In the case of negative Tol
!           the calculations require at least one more iteration (but potentially
!           can result in many more iterations, up to the iteration limit, MaxIter)
!           and the total computation time may increase substantially.
! MaxIter - The maximum number of iterations allowed;
!SpecifNorm - The parameter that determines how the eivenvector is 
!           normalized upon exit. If SpecifNorm=0 then it is normalized 
!           in such a way that xH*B*x=1. If SpecifNorm=1 then xH*x=1. 
!           Otherwise the eigenvector is normalized in such a way that
!           the magnitude of the largest either real or imaginary part
!           is one; 
!  Output parameters:
!       M - A two-dinensional array whose upper triangle (excluding
!           the diagonal) contains the updated matrix LH (of
!           size n) and the lower triangle (including the 
!           diagonal) contains the lower triangle of matrix
!           A-apprlambda*B (of size n); 
!    invD - An array containing the inversed values of
!           the diagonal matrix D (of size n);
!  lambda - The eigenvalue found;
!       x - The eigenvector found;
!  RelAcc - A rough estimate of the relative accuracy reached.
!           It may be somewhat inaccurate. A more accurate estimate
!           would require additional calculations and for the
!           sake of efficiency was not implemented.
! NumIter - The number of iterations performed;
!ErrorCode - An error flag. If ErrorCode=0 then the routine
!           finished succesfully. There may be the following
!           errors on exit:
!           ErrorCode=1 - matrix A-apprlambda*B is singular 
!           or almost singular;
!           ErrorCode=2 - the process did not converge to the 
!           accuracy Tol after the maximum number of 
!           iterations was performed;           

!Arguments
integer         k,n,nM,nB,MaxIter,SpecifNorm,NumIter,ErrorCode
complex(dprec)  M(nM,n),invD(n),B(nB,n),v(n),w(n),x(n)
real(dprec)     apprlambda,Tol,lambda,RelAcc
!Local variables
real(dprec)     NormOfDiff,NormOfDiffPrev,t1,t2,sqrtn
complex(dprec)  tc
logical         notconverged

ErrorCode=0
RelAcc=huge(RelAcc)
NormOfDiffPrev=huge(NormOfDiffPrev)
NumIter=0
!Updating M=L*D*LH factorization up to size n using routine LDLHF
call LDLHF(k,n,M,nM,invD,w,ErrorCode)
if (ErrorCode>0) return
sqrtn=sqrt(2*n*ONE)
!Do inverse iterations until the process converges
!with relative accuracy Tol, or until the number
!of iterations exceeds the limit
notconverged=.true.
do while ((notconverged).AND.(NumIter<MaxIter))
  if (NumIter/=0) v(1:n)=x(1:n)
  call MHMV(n,B,nB,v,w,x)
  call LDLHS(n,M,nM,invD,w,x)
  t1=CMaxAbsReOrIm(n,x)
  call CVScale(n,x,cmplx(1/t1,ZERO,dprec))  
  !tc=(x^{H}v)/(v^{H}v)
  tc=CDotProdQuotient(n,x,v)
  !NormOfDiff=sqrt(2*n)||x-tc*v||
  NormOfDiff=sqrtn*CVDiffEucNorm(n,x,tc,v)
  if (Tol>ZERO) then
    if (NormOfDiff<=Tol) notconverged=.false.
  else
    if ((NormOfDiff>NormOfDiffPrev).AND.(NormOfDiff<=abs(Tol))) notconverged=.false.
    NormOfDiffPrev=NormOfDiff
  endif
  NumIter=NumIter+1
enddo
RelAcc=NormOfDiff
if (notconverged) ErrorCode=2
!Compute x^{H}Bx
t1=VMMHMV(n,B,nB,x)
!Compute x^{H}Mx
t2=VMMHMV(n,M,nM,x)
!Rayleigh quotient
lambda=t2/t1+apprlambda  
select case	(SpecifNorm)
  case (0) !x^{H}Bx=1
    call CVScale(n,x,cmplx(1/sqrt(t1),ZERO,dprec))
  case (1) !x^{H}x=1 
    t1=CDotProdItself(n,x)
    call CVScale(n,x,cmplx(1/sqrt(t1),ZERO,dprec))
endselect

end subroutine GHEPIIS


end module linalg
