F90 = mpiifort
F77 = mpiifort

OPT=-c -O3 -ip -ftz -align all -qopt-report0 -qopt-report-phase=vec
OPT1=-O3 -ip -ftz -align all -qopt-report0 -qopt-report-phase=vec
OPT3=-c -O3 -ip -ftz -align all -qopt-report0 -qopt-report-phase=vec
OPT4=-c -O3 -ip -ftz -align all -qopt-report0 -qopt-report-phase=vec

FILES = globvars.f90 misc.f90 dmng.f linalg.f90 matelem.f90 matform.f90 workproc.f90 main.f90
OBJS = globvars.o misc.o dmng.o linalg.o matelem.o matform.o workproc.o main.o
LIBS = -lblas -llapack

main	 : $(OBJS)
	$(F90) -o main $(OPTL) $(OBJS) $(LIBS) 

main.o	 : $(FILES)
	$(F90) -c $(OPT) $?

globvars.o : globvars.f90
	$(F90) -c $(OPT) $(LIBS) $?

misc.o : misc.f90
	$(F90) -c $(OPT) $?

dmng.o : dmng.f
	$(F77) -c $(OPT) $?

linalg.o : linalg.f90 
	$(F90) -c $(OPT) $?

matelem.o : matelem.f90 globvars.f90
	$(F90) -c $(OPT) $?

matform.o : matform.f90 matelem.f90 
	$(F90) -c $(OPT) $?

workproc.o : workproc.f90 
	$(F90) -c $(OPT) $?
