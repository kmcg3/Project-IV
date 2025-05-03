MODULE potential
  use global_variables
  implicit none

  ! Functions to override
  logical, parameter :: run_init = .true.
  logical, parameter :: perturb_coordinates_override = .true.

  ! Model global variables
  DOUBLE PRECISION, ALLOCATABLE :: COORDS(:,:), K_ARRAY(:), R_ARRAY(:), B_ARRAY(:)
  DOUBLE PRECISION, ALLOCATABLE :: COORDS_0(:,:), ANGLES(:), LENGTHS(:), GRAD(:), GRAD_A(:), GRAD_L(:), GRAD_SPAT(:,:)
  INTEGER, ALLOCATABLE :: BONDS(:,:), TRIANGLES(:,:)
  DOUBLE PRECISION :: F_I, F_S, L, LS_0, H, K1, K2, K3, R1, R2, R3, BI, BE, HS, CM1, KC, NORM_ANG, KY, WY,LOWER_RAT, KL
  DOUBLE PRECISION :: KAREA, CAREA, AREA0
  INTEGER :: N_ROD, N_SEG, i, j, k, N_BONDS, IO

!LJ STRETCHING, LENGTHS FIXED NOT THE SAME, VARYING BENDING ENERGY ALONG THE LENGTH (B_ARRAY(:,:)),

CONTAINS
!------------------------------------------------------------------------------------------
  !Wrapper subroutines
  subroutine init()
    ! Wrapped for code that needs to be called
    write(*,*) "InitWetting"
    call INIT_SYSTEM()
    write(*,*) "InitWetting end"
  end subroutine

  subroutine calc_energy_gradient()
    implicit none
    ! Wrapper to phase field model
    call IMPLEMENT_POTENTIAL(X, G, E, .true.)
   
  end subroutine
!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
!System initialisation
SUBROUTINE INIT_SYSTEM()
  IMPLICIT NONE
    INTEGER :: CUR, NODE1, NODE2, NODE0, S1

  DOUBLE PRECISION :: PI, SUM_TEMP
  character(len=10)       :: datechar,timechar,zonechar
  integer                 :: values(8),itime1
  CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
  itime1= values(7)*39 + values(8)
  CALL SDPRND(itime1)

 !Define the og segment lengths and vertical separation 
  LS_0 = L/(N_SEG)
  HS = H/(N_ROD-1)
  
  ALLOCATE(COORDS(N_ROD*(N_SEG+1),2), K_ARRAY(3), R_ARRAY(3), COORDS_0(N_ROD,2), ANGLES(N_ROD*(N_SEG+1)), LENGTHS(N_ROD*(N_SEG+1)), GRAD(N), GRAD_A(N/2), GRAD_L(N/2), GRAD_SPAT(N_ROD*(N_SEG+1),2))
  COORDS=0.0
  ANGLES=0.0
  LENGTHS=LS_0 !???!
  GRAD=0.0
  COORDS_0=0.0




  !Fill the elastic constant and equilibrium bond length arrays
  K_ARRAY(1)=K1
  K_ARRAY(2)=K2
  K_ARRAY(3)=K3
  R_ARRAY(1)=R1
  R_ARRAY(2)=R2
  R_ARRAY(3)=R3

  !Read in the coordinates of the rod starting points
  OPEN(1,FILE='initpos')
  DO i=1,N_ROD
  	READ(1,*) COORDS_0(i,:)
  ENDDO
  CLOSE(1)

  !Read in the bonds
  N_BONDS=0
  OPEN(1, FILE='bonds')
  DO
  READ(1,*, iostat=io)
  IF  (io/=0) EXIT
  	N_BONDS=N_BONDS+1
  END DO
  CLOSE(1)
  
  ALLOCATE(BONDS(N_BONDS,3))

  OPEN(1,FILE='bonds')
  DO i=1,N_BONDS
  	READ(1,*) BONDS(i,:)
  ENDDO
  CLOSE(1)

!Make a list of the node triplets used in the area repulsion term
  AREA0=0.5*LS_0*HS !Equilibrium area of triangles (for now they are uniform)
  ALLOCATE(TRIANGLES(2*(N_ROD-1)*N_SEG,3))
  TRIANGLES=0
  S1=1
 ! DO i=1,N_ROD-1
!	DO j=1,N_SEG+1
!		IF (j<N_SEG+1) THEN
!			NODE0=(i-1)*(N_SEG+1)+j
!			NODE1=(i)*(N_SEG+1)+j+1
!			NODE2=(i)*(N_SEG+1)+j
!			TRIANGLES(S1,1)=NODE0
!			TRIANGLES(S1,2)=NODE1
!			TRIANGLES(S1,3)=NODE2	
!			S1=S1+1
!		ENDIF
!		
!		IF (j>1) THEN
!			NODE0=(i-1)*(N_SEG+1)+j
!			NODE1=(i)*(N_SEG+1)+j
!			NODE2=(i)*(N_SEG+1)+j-1
!			TRIANGLES(S1,1)=NODE0
!			TRIANGLES(S1,2)=NODE1
!			TRIANGLES(S1,3)=NODE2	
!			S1=S1+1
!		ENDIF
!		
!	ENDDO
!  ENDDO

  DO i=1,N_ROD-1
	DO j=1,N_SEG
		NODE0=(i-1)*(N_SEG+1)+j
		NODE1=(i)*(N_SEG+1)+j+1
		NODE2=(i)*(N_SEG+1)+j
		TRIANGLES(S1,1)=NODE0
		TRIANGLES(S1,2)=NODE1
		TRIANGLES(S1,3)=NODE2
		S1=S1+1	
	ENDDO
  ENDDO
			
  DO i=2,N_ROD
	DO j=2,N_SEG+1
		NODE0=(i-1)*(N_SEG+1)+j
		NODE1=(i-2)*(N_SEG+1)+j-1
		NODE2=(i-2)*(N_SEG+1)+j
		TRIANGLES(S1,1)=NODE0
		TRIANGLES(S1,2)=NODE1
		TRIANGLES(S1,3)=NODE2
		S1=S1+1	
	ENDDO
  ENDDO

!DO i=1,2*(N_ROD-1)*N_SEG
!	PRINT *, TRIANGLES(i,:)
!ENDDO
!STOP


END SUBROUTINE INIT_SYSTEM
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!System perturbation
SUBROUTINE PERTURB_COORDINATES() !THIS DOES NOTHING?
  IMPLICIT NONE
  INTEGER :: i, j, CUR 

  !This is where we can make initial changes to the angles from the initial input angles X
  ANGLES=X(1:(N/2)) !FINE
  LENGTHS=X((N/2)+1:N) !FINE
  !PRINT *, "perturb_coords, LENGTHS = " ,LENGTHS
 ! PRINT *, "X = ", X((N/2)+1:N)
!PRINT *, "LS0=", LS_0
!PRINT *, "LENGTHS=", LENGTHS
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		IF (j <= N_SEG-25) THEN
			!ANGLES(CUR) = 0.00
		ENDIF
	ENDDO
  ENDDO

  X(1:(N/2))=ANGLES !FINE
  X((N/2)+1:N)=LENGTHS
  !PRINT *, "After PERTURB_COORDS X = " , X
!   !PRINT *, "X = ", X((N/2)+1:N)
! OPEN(1,FILE='coords.perturb')
! WRITE(1,'(F20.10)') X
! CLOSE(1)

! OPEN(1,FILE='coords.perturb')
! WRITE(1,'(F20.10)') X !DOES THIS WORK??? 
! CLOSE(1)



END SUBROUTINE PERTURB_COORDINATES

!NEW...
! SUBROUTINE PERTURB_LENGTHS()
!   IMPLICIT NONE
!   INTEGER :: i, j, CUR 

!   !This is where we can make initial changes to the LENGTHS from the initial input LENGTHS X
!   !X((N/2)+1:N)=0.0001
!   LENGTHS=X((N/2)+1:N)
!   PRINT *, "PERT COORDS LENGTHS = " , LENGTHS

!   DO i=1,N_ROD
! 	DO j=1,N_SEG-1
! 		CUR=(i-1)*(N_SEG-1)+j
! 		IF (j <= N_SEG-25) THEN
! 			!LENGTHS(CUR) = 0.00
! 		ENDIF
! 	ENDDO
!   ENDDO

!   X((N/2)+1:N)=LENGTHS

! OPEN(1,FILE='coords.perturb', STATUS='UNKNOWN', POSITION='APPEND') !CHANGED THIS...NOT SURE IT WILL WORK
! WRITE(1,'(F20.10)') LENGTHS
! CLOSE(1)



!END SUBROUTINE PERTURB_LENGTHS


!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!Compute the energy, energy derivatives and gradient
Subroutine IMPLEMENT_POTENTIAL(X,G,E,GTEST) 
  IMPLICIT NONE
  LOGICAL GTEST
  DOUBLE PRECISION :: X(N), E, G(N)
!PRINT *, "X", X
ANGLES=X(1:(N/2))
LENGTHS=X((N/2)+1:N)

  CALL UNPACK_COORDINATES()

  E = 0.0
  GRAD_A = 0.0
  GRAD_L = 0.0

  CALL COMPUTE_ENERGY(E)

  IF (GTEST) THEN
   	CALL COMPUTE_GRAD()
  ENDIF


  G(1:N/2) = GRAD_A * CM1
  G((N/2)+1:N) = GRAD_L * CM1

 ! G = GRAD * CM1
  E = E * CM1
  X(1:(N/2))=ANGLES !FINE
  X((N/2)+1:N)=LENGTHS
  !PRINT *, "After PERTURB_COORDS X = " , X
  !PRINT *, "X = ", X((N/2)+1:N)
! OPEN(1,FILE='coords.perturb')
! WRITE(1,'(F20.10)') X
! CLOSE(1)
! CALL TEST_GRADIENT(G,X)

END SUBROUTINE IMPLEMENT_POTENTIAL
!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
!Unpack the angles to form the vertex coordinates  
SUBROUTINE UNPACK_COORDINATES()
  IMPLICIT NONE
  INTEGER :: CUR
  !DOUBLE PRECISION :: LENGTHS(N/2), ANGLES(N/2)
  !Coordinates of 1st and 2nd vertex of each rod
       !PRINT *, "In UNPACK, LENGTHS = , LENGTHS
  !PRINT *, "IN UNPACK ANGLES=", ANGLES


  DO i=1,N_ROD

    COORDS((i-1)*(N_SEG+1)+1,1) = COORDS_0(i,1);
    COORDS((i-1)*(N_SEG+1)+1,2) = COORDS_0(i,2); 
    
    COORDS((i-1)*(N_SEG+1)+2,1) = COORDS((i-1)*(N_SEG+1)+1,1) + LENGTHS((i-1)*N_SEG+1); !check this lengths 
    COORDS((i-1)*(N_SEG+1)+2,2) = COORDS((i-1)*(N_SEG+1)+1,2);  
  ENDDO

  !Coordinates of remaining vertices of rod 1 to (N_ROD-1)
  DO i=1,N_ROD
  DO j=3,N_SEG+1
		CUR=(i-1)*(N_SEG-1)+j
        COORDS((i-1)*(N_SEG+1)+j,1) = COORDS((i-1)*(N_SEG+1)+j-1,1) + LENGTHS(CUR-2)*cos(ANGLES((i-1)*(N_SEG-1) + j-2)); 
        COORDS((i-1)*(N_SEG+1)+j,2) = COORDS((i-1)*(N_SEG+1)+j-1,2) + LENGTHS(CUR-2)*sin(ANGLES((i-1)*(N_SEG-1) + j-2));
  ENDDO
  ENDDO
END SUBROUTINE UNPACK_COORDINATES
!------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------
!Energy computation
SUBROUTINE COMPUTE_ENERGY(E)
  IMPLICIT NONE
  DOUBLE PRECISION :: E, EB, ES, EY, ELI, ELS, ETH, EC, EAREA, ETEMP, EL, COORD1(2), COORD2(2), R, NORMAL(2), TEMP
  DOUBLE PRECISION :: DIFF_0, COS_DIFF, VERTDIF,X0,XP,XM,Y0,YP,YM,AREA 
  INTEGER :: CUR, NODE1, NODE2, BONDTYPE, CURP, CURM, CUR2P, CUR2M, CUR2

  E=0.0
  EB=0.0
  ES=0.0
  ELI=0.0
  ELS=0.0
  EC=0.0
  EY=0.0
  ETH=0.0
  EL=0.0


!PRINT *, "X ENG", X(1)
!PRINT *, "lengths and angles in energy = ", LENGTHS, ANGLES
!PRINT *, LENGTHS
! DO i=1,N_ROD
!	DO j=1,N_SEG+1
!		CUR=(i-1)*(N_SEG+1)+j
!		PRINT *, COORDS(CUR,:)
!	ENDDO
! ENDDO
 !PRINT *, "X Energy", X(1)
 !PRINT *, "L Energy", LENGTHS(1)
 ! PRINT *, "A Energy", ANGLES(1)

  !1) Bending energy
  !1.1) Interior rods
  DO i=2,N_ROD-1
	DO j=2,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
    
    EB = EB + (B_ARRAY(i)/LENGTHS(CUR))*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2

! !NEW CURVATURE -------------
!   D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!   ADIFF = ANGLES(CUR)-ANGLES(CUR-1)
  
!   EB = EB + B_ARRAY(i)*8*LENGTHS(CUR)*(sin(ADIFF)**2)/D
!   !-------------------------
	ENDDO
  ENDDO
  EB = EB

  !1.2) Exterior rods
  DO j=2,N_SEG-1
	CUR=j
	EB = EB + (2*BE/LENGTHS(CUR))*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2

  !NEW CURVATURE -------------
!   D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!   ADIFF = ANGLES(CUR)-ANGLES(CUR-1)
  
!   EB = EB + BE*8*LENGTHS(CUR)*(sin(ADIFF)**2)/D
!   !-------------------------

  ENDDO
 ! PRINT *, "EB", EB

  DO j=2,N_SEG-1
	CUR=(N_ROD-1)*(N_SEG-1)+j
	EB = EB +(2*BE/LENGTHS(CUR))*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2
! !NEW CURVATURE -------------
!   D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!   ADIFF = ANGLES(CUR)-ANGLES(CUR-1)
  
!   EB = EB + BE*8*LENGTHS(CUR)*(sin(ADIFF)**2)/D
!   !-------------------------

  ENDDO
  
  !2) Stretching energy
  DO i=1,N_BONDS
	NODE1=BONDS(i,1)
	NODE2=BONDS(i,2)
	BONDTYPE=BONDS(i,3)

	R=SQRT( (COORDS(NODE1,1)-COORDS(NODE2,1))**2 + (COORDS(NODE1,2)-COORDS(NODE2,2))**2 )

	!With 1/R term
	!ES = ES + K_ARRAY(BONDTYPE)*(R-R_ARRAY(BONDTYPE))**2*(1.0/R)
	!Without 1/R term
	!ES = ES + K_ARRAY(BONDTYPE)*(R-R_ARRAY(BONDTYPE))**2 OG

  !Harmonic ------------
  ES = ES + K_ARRAY(BONDTYPE)*((R-R_ARRAY(BONDTYPE))/R_ARRAY(BONDTYPE))**2 !

  !LJ ------------------
  !ES = ES + K_ARRAY(BONDTYPE)*((R_ARRAY(BONDTYPE)/R)**12-(R_ARRAY(BONDTYPE)/R)**6) !K^2? (EPSILON) AND *LENGTHS???
  ENDDO
  ES = ES !DO I WANT THIS STILL?
! PRINT *, "ES", ES
  
  !3) Interocular load
  DO i=1,N_ROD
	CUR=i*(N_SEG+1)
	ELI = ELI + COORDS(CUR,1)
  ENDDO
  ELI=ELI*F_I
!PRINT *, "ELI", ELI

  !4) Sclera load

  ETEMP=0.0
  DO j=2,N_SEG+1
	CUR=j
	ETEMP = ETEMP + COORDS(CUR,2)-COORDS(1,2)
  !PRINT *, "ETEMP", ETEMP
  ENDDO
  ELS = ELS-ETEMP*F_S
  ! PRINT *, "ETE", ETEMP*F_S
  ! PRINT *, "ELS", ELS
  !ELS is just 0 each time

  ETEMP=0.0
  DO j=2,N_SEG+1
  	CUR=(N_ROD-1)*(N_SEG+1)+j
  	ETEMP = ETEMP + COORDS(CUR,2)-COORDS((N_ROD-1)*(N_SEG+1)+1,2)
   ! PRINT *, "C", COORDS(CUR,2)-COORDS((N_ROD-1)*(N_SEG+1)+1,2)
   ! PRINT *, "E", ETEMP
  ENDDO
  !  PRINT *, "ete", ETEMP
  !  PRINT *, "ETEMP", ETEMP
  ELS = ELS+ETEMP*F_S
  ! PRINT *, "F_S", F_S
  ! PRINT *, "ELS", ELS

  !5) Angular soft constraints of rod ends
  !5.1) TWO ENDS SHOULD BE ALIGNED
  DO i=1,N_ROD
  CUR=(N_ROD-1)*(N_SEG-1)+j
  VERTDIF =COORDS((i-1)*(N_SEG+1)+2,2)-COORDS((i-1)*(N_SEG+1)+N_SEG+1,2)
  EC = EC + 0.5*KC*VERTDIF*VERTDIF
  
  ENDDO

 !PRINT *, "EC", EC
  	  !PRINT *, "lengths in eng3 = ", LENGTHS(i)
  !5.2) END ANGLES CONSTRAINED (PINNED)
  DO i=1,N_ROD
    CUR=(i-1)*(N_SEG-1)+1
    EC = EC + 0.5*KC*(1-COS(ANGLES(CUR)))*(1-COS(ANGLES(CUR)))
    ENDDO
    
    DO i=1,N_ROD
    CUR=(i-1)*(N_SEG-1)+N_SEG-1
    EC = EC + 0.5*KC*(1-COS(ANGLES(CUR)))*(1-COS(ANGLES(CUR)))
    ENDDO
  ! PRINT *, "EC", EC
  !6) Energy penalising a node on a rod below becoming above the same node on the rod above
  ! DO i=1,N_ROD-1
	! DO j=1,N_SEG+1
	! 	CUR = (i-1)*(N_SEG+1) + j
	! 	CURP = (i)*(N_SEG+1) + j
	! 	TEMP= 0.0
	! 	TEMP = ( 1.0 - TANH((COORDS(CURP,2)-COORDS(CUR,2))/WY) )
	! 	TEMP = TEMP*( 2.0 - (COORDS(CURP,2)-COORDS(CUR,2)) )
	! 	TEMP = TEMP * KY * 0.5
		
	! EY = EY+TEMP

	! ENDDO
  ! ENDDO OG

   !6) Repulsion term that prevents the node triangles having <0 area
  DO i=1,2*(N_ROD-1)*N_SEG
	X0=COORDS(TRIANGLES(i,1),1)
	XP=COORDS(TRIANGLES(i,2),1)
	XM=COORDS(TRIANGLES(i,3),1)
	Y0=COORDS(TRIANGLES(i,1),2)
	YP=COORDS(TRIANGLES(i,2),2)
	YM=COORDS(TRIANGLES(i,3),2)
	AREA=0.5*((XP-X0)*(YM-Y0)-(XM-X0)*(YP-Y0))
!PRINT *, "AREA", AREA
!PRINT *, "AREA0", AREA0
	EAREA=EAREA+KAREA*EXP(-CAREA*AREA/AREA0)
  ENDDO
  !PRINT *, "EAREA", EAREA
 
  !7) Energy penalising an angle geeing above pi
  
  !DO i=1,N_ROD
	!DO j=1,N_SEG-1
	!	CUR=(i-1)*(N_SEG-1)+j
  ! TEMP = ANGLES(CUR)-3.14159265359
	!	ETH = ETH + KY*TANH((ABS(ANGLES(CUR))-3.14159265359)/0.3);
	!ENDDO
  !ENDDO
  
!   !8) Energy penalising segment lengths being zero
!HARMONIC--------------
  DO i=1,N_ROD
  DO j=1,N_SEG+1
    CUR = (i-1)*(N_SEG-1) + j
!     IF (LENGTHS(CUR) <= 0) THEN
!    PRINT *, "LENGTHS is negative (or zero) ", CUR
! END IF
    EL = EL + (KL/2)*((LENGTHS(CUR)-LS_0)**2)
  ENDDO
  ENDDO

! !LJ --------------------
!   DO i=1,N_ROD
!   DO j=1,N_SEG+1
!     CUR = (i-1)*(N_SEG-1) + j
!     EL = EL + (KL**2)*((LS_0/LENGTHS(CUR))**12-(LS_0/LENGTHS(CUR))**6)
!   ENDDO
!   ENDDO
!PRINT *, "E", EL
! !  E = EB + ES + ELI + ELS + EC +EY

E = EB + ES + ELI + EC + EAREA + EL

!PRINT *, EB, ES, ELI, ELS, EC, EY, ETH, EL, 
!PRINT *, EAREA
!PRINT *, "LENGTHS 1", LENGTHS(1)
!PRINT *,"E", E
!PRINT *, "Energy computed"

END SUBROUTINE COMPUTE_ENERGY
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!Gradient computation
SUBROUTINE COMPUTE_GRAD()
  IMPLICIT NONE
  DOUBLE PRECISION :: THETAM, THETAP, TEMP, R, WEIGHTX, WEIGHTY, LENGTHSM, LENGTHSP, NUMM, NUMP, NUMML, NUMPL
  DOUBLE PRECISION :: DIFF_0, DIFF_P, DIFF_M, DIFF_0_SQ, DIFF_P_SQ, DIFF_M_SQ, DIFF_0_Y, DIFF_P_Y, DIFF_M_Y
  DOUBLE PRECISION :: TEMP_P, TEMP_M, COS_DIFF, GTEMP
  DOUBLE PRECISION :: X0, XP, XM, Y0, YP, YM, AREA, VPRE
  INTEGER :: CUR, CURP, CURM,  NODE1, NODE2, BONDTYPE, NODE_J, NODE_K, END_NODE_K, CUR2, CUR2P, CUR2M, CUR2PP, CUR2MM, GSIGN

  GRAD_A=0.0
  GRAD_L=0.0
  GRAD_SPAT=0.0


 !1) Bending gradient 
  !1.1) Internal rods, not including ends
  !1.1.1) GRAD_A
  DO i=2,N_ROD-1
	DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

			GRAD_A(CUR) = GRAD_A(CUR) + 2*B_ARRAY(i)/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2))-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  

ENDDO 
ENDDO
  
  !1.1.2) GRAD_L
DO i=2,N_ROD-1
	DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

    THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))


		GRAD_L(CUR) = GRAD_L(CUR) - 2*B_ARRAY(i)/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CUR-1))) )**2

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D)+(NUMPL)/(D))
!  !-------------------------  

ENDDO 
ENDDO


  !1.2) Internal rods, first segment (angle 1 = 0)  
  !1.2.1) GRAD_A
  DO i=2,N_ROD-1
	j=1
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))
 
		GRAD_A(CUR) = GRAD_A(CUR) + 2*B_ARRAY(i)/LENGTHS(CUR)*( - TAN(THETAP)/(COS(THETAP))**2 )

!     ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/D)*(-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  
  ENDDO

  !1.2.2) GRAD_L
  DO i=2,N_ROD-1
    j=1
      CUR = (i-1)*(N_SEG-1) + j
      CURP = (i-1)*(N_SEG-1) + j+1
  
      THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))
  
      GRAD_L(CUR) = GRAD_L(CUR) - 2*B_ARRAY(i)/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CURP)-ANGLES(CUR))) )**2

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/(D**2))*((NUMPL)/(D))
!  !-------------------------  

  ENDDO


  !1.3) Internal rods, last segment (no N_SEG+1 contribution)
  !1.3.1) GRAD_A
  DO i=2,N_ROD-1
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
 
		GRAD_A(CUR) = GRAD_A(CUR) + 2*B_ARRAY(i)/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 )

!     ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2)))
!  !-------------------------  
ENDDO

	!1.3.2) GRAD_L
DO i=2,N_ROD-1
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))

		GRAD_L(CUR) = GRAD_L(CUR) - 2*B_ARRAY(i)/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CURM))) )**2

    ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*B_ARRAY(i)*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D)+)
!  !-------------------------  


ENDDO

  !1.4) External rods, not including ends
  !1.4.1) GRAD_A
  i=1
DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )

!     ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2))-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  

ENDDO 
 
  i= N_ROD
DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 - TAN(THETAP)/(COS(THETAP))**2 )

!     ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2))-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  
ENDDO 



!1.4.2) GRAD_L
i=1
DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

	GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CURM))) )**2

    ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D)+(NUMPL)/(D))
!  !-------------------------  

ENDDO 
 
  i= N_ROD
DO j=2,N_SEG-2
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))
		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

  GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CURM))) )**2 

    ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D)+(NUMPL)/(D))
!  !-------------------------  

ENDDO 

  !1.5) External rods, first segments (angle 1 = 0)
  i=1
	j=1
  
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))

		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( - TAN(THETAP)/(COS(THETAP))**2 )

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  
  GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CURP)-ANGLES(CUR))) )**2

   ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*((NUMPL)/(D))
!  !-------------------------  
 


  i=1
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))


		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 ) 

!     ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2)))
!  !-------------------------  

    GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CURM))) )**2

   ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D))
!  !-------------------------  

  
   !1.4) External rods,  last segment (no N_SEG+1 contribution)

   i=N_ROD
	j=1
		CUR = (i-1)*(N_SEG-1) + j
		CURP = (i-1)*(N_SEG-1) + j+1

		THETAP=0.5*(ANGLES(CURP)-ANGLES(CUR))


		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( - TAN(THETAP)/(COS(THETAP))**2 )

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(-SIN(4*THETAP)-(NUMP/(D**2)))
!  !-------------------------  

  GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CURP)-ANGLES(CUR))) )**2

    ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*((NUMPL)/(D))
!  !-------------------------  



  !1.4) External rods,  last segment (no N_SEG+1 contribution)
  i=N_ROD
	j=N_SEG-1
		CUR = (i-1)*(N_SEG-1) + j
		CURM = (i-1)*(N_SEG-1) + j-1

		THETAM=0.5*(ANGLES(CUR)-ANGLES(CURM))

		GRAD_A(CUR) = GRAD_A(CUR) + 2*BE/LENGTHS(CUR)*( TAN(THETAM)/(COS(THETAM))**2 ) 

! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMM = (SIN(2*THETAM)**3)*LENGTHS(CUR)*LENGTHS(CURM)
!     NUMP = (SIN(2*THETAP)**3)*LENGTHS(CUR)*LENGTHS(CURP)

!     GRAD_A(CUR) = GRAD_A(CUR) + 8*BE*LENGTHS(CUR)*(1/D)*(SIN(4*THETAM)+(NUMM/(D**2)))
!  !-------------------------  

  GRAD_L(CUR) = GRAD_L(CUR) - 2*BE/(LENGTHS(CUR)**2)*( tan(0.5*(ANGLES(CUR)-ANGLES(CURM))) )**2

    ! ! !NEW CURVATURE -------------
!     D = ((COORDS(CUR+1,1)-COORDS(CUR-1,1))**2 + (COORDS(CUR+1,2)-COORDS(CUR-1,2))**2)**(1/2)
!     NUMML = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURM)*COS(2*THETAM))
!     NUMPL = LENGTHS(CUR)*(LENGTHS(CUR)+LENGTHS(CURP)*COS(2*THETAP))

!     GRAD_L(CUR) = GRAD_L(CUR) + 8*BE*LENGTHS(CUR)*(1/(D**2))*(-(NUMML)/(D))
!  !-------------------------  


!   !2) Stretching gradient
  !2.1) Start by computing the dE/dxi and dE/dyi
  DO i=1,N_BONDS
	NODE1=BONDS(i,1)
	NODE2=BONDS(i,2)
	BONDTYPE=BONDS(i,3)

	R=SQRT( (COORDS(NODE1,1)-COORDS(NODE2,1))**2 + (COORDS(NODE1,2)-COORDS(NODE2,2))**2 )
	
	! X-derivatives
	!With 1/R term
	!TEMP=K_ARRAY(BONDTYPE)*(COORDS(NODE1,1)-COORDS(NODE2,1))*(1.0-R_ARRAY(BONDTYPE)/R)*(1/R*(1-0.5*(1.0-R_ARRAY(BONDTYPE)/R)))
	!Without 1/R term
	!TEMP=K_ARRAY(BONDTYPE)*(COORDS(NODE1,1)-COORDS(NODE2,1))*(1.0-R_ARRAY(BONDTYPE)/R) OG

  !HARMONIC -------------
  TEMP=2*(K_ARRAY(BONDTYPE)/(R_ARRAY(BONDTYPE)**2))*(COORDS(NODE1,1)-COORDS(NODE2,1))*(1.0-R_ARRAY(BONDTYPE)/R) !PROPER
    !--------------------

  !LJ -------------------
  !TEMP=(K_ARRAY(BONDTYPE)*(((-12/(R**13))*((R_ARRAY(BONDTYPE))**12)+((6/(R**7))*((R_ARRAY(BONDTYPE)**6))))*((COORDS(NODE1,1)-COORDS(NODE2,1))/R)))
    !--------------------
	GRAD_SPAT(NODE1,1) = GRAD_SPAT(NODE1,1) + TEMP
	GRAD_SPAT(NODE2,1) = GRAD_SPAT(NODE2,1) - TEMP



	! Y-derivatives
	!With 1/R term
	!TEMP=K_ARRAY(BONDTYPE)*(COORDS(NODE1,2)-COORDS(NODE2,2))*(1.0-R_ARRAY(BONDTYPE)/R)*(1/R*(1-0.5*(1.0-R_ARRAY(BONDTYPE)/R)))
	!Without 1/R term
	!TEMP=K_ARRAY(BONDTYPE)*(COORDS(NODE1,2)-COORDS(NODE2,2))*(1.0-R_ARRAY(BONDTYPE)/R) OG

  !HARMONIC -----------
  TEMP=2*(K_ARRAY(BONDTYPE)/(R_ARRAY(BONDTYPE)**2))*(COORDS(NODE1,2)-COORDS(NODE2,2))*(1.0-R_ARRAY(BONDTYPE)/R) !PROPER
  !--------------------

  !LJ -----------------
  !TEMP=(K_ARRAY(BONDTYPE)*(((-12/(R**13))*((R_ARRAY(BONDTYPE))**12)+((6/(R**7))*((R_ARRAY(BONDTYPE)**6))))*((COORDS(NODE1,2)-COORDS(NODE2,2))/R)))
  !--------------------
	GRAD_SPAT(NODE1,2) = GRAD_SPAT(NODE1,2) + TEMP
	GRAD_SPAT(NODE2,2) = GRAD_SPAT(NODE2,2) - TEMP
  ENDDO


  !2.2) Now convert these to dE/dtheta_i
  !2.2.1) Interior rods 
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		! 1D location of current angle
		CUR=(i-1)*(N_SEG-1)+j
	
	DO k=j+1,N_SEG+1
			! 1D location of a node which the current angle affects
			NODE_K=(i-1)*(N_SEG+1)+k+1
			
			GRAD_A(CUR) = GRAD_A(CUR) - LENGTHS(CUR)*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LENGTHS(CUR)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
	ENDDO		
	ENDDO
  ENDDO

  !2.3) dE/dl_i
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		! 1D location of current length
		CUR=(i-1)*(N_SEG-1)+j-1
	
	DO k=j+1,N_SEG+1
			! 1D location of a node which the current length affects
			NODE_K=(i-1)*(N_SEG+1)+k !is this right with the +1?
			
			GRAD_L(CUR) = GRAD_L(CUR) + COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
	ENDDO		
	ENDDO
  ENDDO


  !3) Interocular load
  !3.1)Internal rods
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		CUR=(i-1)*(N_SEG-1)+j
		GRAD_A(CUR) = GRAD_A(CUR) - F_I*LENGTHS(CUR)*SIN(ANGLES(CUR)) 
    GRAD_L(CUR) = GRAD_L(CUR) +F_I*COS(ANGLES(CUR))
    !PRINT *,'GRAD', GRAD_A(CUR)
   !PRINT *, 'F_I', F_I
	ENDDO
  ENDDO

  ! !4) Sclera load
  ! i=1
	! DO j=1,N_SEG-1
	! 	CUR=(i-1)*(N_SEG-1)+j
	! 	GRAD_A(CUR) = GRAD_A(CUR) - F_S*LENGTHS(CUR)*COS(ANGLES(CUR))*(N_SEG-j)
  !   !PRINT *, "FS", F_S
  !   GRAD_L(CUR) = GRAD_L(CUR) - F_S*SIN(ANGLES(CUR))*(N_SEG-j)
	! ENDDO
 
  !  i=N_ROD
	! DO j=1,N_SEG-1
	! 	CUR=(i-1)*(N_SEG-1)+j
	! 	GRAD_A(CUR) = GRAD_A(CUR) + F_S*LENGTHS(CUR)*COS(ANGLES(CUR))*(N_SEG-j)
  !  GRAD_L(CUR) = GRAD_L(CUR) + F_S*SIN(ANGLES(CUR))*(N_SEG-j)
	! ENDDO

!5)Angular constraint
 ! 5.1) End aligment
  DO i=1,N_ROD
    DO j=1,N_SEG-1
      CUR=(i-1)*(N_SEG-1)+j 
      GRAD_A(CUR) = GRAD_A(CUR) - KC*(COORDS((i-1)*(N_SEG+1)+2,2)-COORDS((i-1)*(N_SEG+1)+N_SEG+1,2))*LENGTHS(CUR)*COS(ANGLES(CUR))
      GRAD_L(CUR) = GRAD_L(CUR) - KC*(COORDS((i-1)*(N_SEG+1)+2,2)-COORDS((i-1)*(N_SEG+1)+N_SEG+1,2))*SIN(ANGLES(CUR))
    ENDDO
  ENDDO
!   !5.2) END ANGLES CONSTRAINED (PINNED)
  DO i=1,N_ROD
    CUR=(i-1)*(N_SEG-1)+1
    GRAD_A(CUR) =  GRAD_A(CUR)+KC*SIN(ANGLES(CUR))*(1-COS(ANGLES(CUR)))
   ENDDO
 
   DO i=1,N_ROD
    CUR=(i-1)*(N_SEG-1)+N_SEG-1
    GRAD_A(CUR) =  GRAD_A(CUR)+KC*SIN(ANGLES(CUR))*(1-COS(ANGLES(CUR)))
   ENDDO	
 
  
!6) Gradient of overlap penalty

  !6) Repulsion term that prevents the node triangles having <0 area
  GRAD_SPAT=0.0
  DO i=1,2*(N_ROD-1)*N_SEG
	X0=COORDS(TRIANGLES(i,1),1)
	XP=COORDS(TRIANGLES(i,2),1)
	XM=COORDS(TRIANGLES(i,3),1)
	Y0=COORDS(TRIANGLES(i,1),2)
	YP=COORDS(TRIANGLES(i,2),2)
	YM=COORDS(TRIANGLES(i,3),2)

	AREA=0.5*((XP-X0)*(YM-Y0)-(XM-X0)*(YP-Y0))

	VPRE=EXP(-CAREA*AREA/AREA0)


	GRAD_SPAT(TRIANGLES(i,1),1)=GRAD_SPAT(TRIANGLES(i,1),1)+VPRE*(-YM+YP)
	GRAD_SPAT(TRIANGLES(i,2),1)=GRAD_SPAT(TRIANGLES(i,2),1)+VPRE*(YM-Y0)
	GRAD_SPAT(TRIANGLES(i,3),1)=GRAD_SPAT(TRIANGLES(i,3),1)+VPRE*(-YP+Y0)

	GRAD_SPAT(TRIANGLES(i,1),2)=GRAD_SPAT(TRIANGLES(i,1),2)-VPRE*(-XM+XP)
	GRAD_SPAT(TRIANGLES(i,2),2)=GRAD_SPAT(TRIANGLES(i,2),2)-VPRE*(XM-X0)
	GRAD_SPAT(TRIANGLES(i,3),2)=GRAD_SPAT(TRIANGLES(i,3),2)-VPRE*(-XP+X0)
  ENDDO
  GRAD_SPAT=GRAD_SPAT*(-CAREA*KAREA/AREA0*0.5)


  ! !6.1) Compute the spatial gradients OG
  ! !6.1.1) Interior rods
  ! DO i=2,N_ROD-1
	! DO j=1,N_SEG+1
	! 	CUR = (i-1)*(N_SEG+1) + j
	! 	CURP = (i)*(N_SEG+1) +j
	! 	CURM= (i-2)*(N_SEG+1) +j
		
	! 	TEMP = 1.0 - tanh( (COORDS(CURP,2)-COORDS(CUR,2))/WY )
	! 	TEMP = TEMP + (2.0-COORDS(CURP,2)+COORDS(CUR,2))*( 1.0/WY * 1.0/(COSH((COORDS(CURP,2)-COORDS(CUR,2))/WY))**2 )
	! 	TEMP = TEMP*KY*0.5

	! 	TEMP_P = -1.0 + tanh( (COORDS(CUR,2)-COORDS(CURM,2))/WY )
	! 	TEMP_P = TEMP_P - (2-COORDS(CUR,2)+COORDS(CURM,2))*( 1.0/WY * 1.0/(COSH((COORDS(CUR,2)-COORDS(CURM,2))/WY))**2 )
	! 	TEMP_P = TEMP_P*KY*0.5

  !  		GRAD_SPAT(CUR,2)=TEMP+TEMP_P

	! ENDDO
  ! ENDDO

  ! !6.1.2) Exterior rods OG
  ! !Rod 1
  ! DO j=1,N_SEG+1
	! CUR = (1-1)*(N_SEG+1) + j
	! CURP = (1)*(N_SEG+1) +j

	! TEMP = 1.0 - tanh( (COORDS(CURP,2)-COORDS(CUR,2))/WY )
	! TEMP = TEMP + (2-COORDS(CURP,2)+COORDS(CUR,2))*( 1.0/WY * 1.0/(COSH((COORDS(CURP,2)-COORDS(CUR,2))/WY))**2 )
	! TEMP = TEMP*KY*0.5

  !  	GRAD_SPAT(CUR,2)=TEMP
  ! ENDDO
  ! !Rod N
  ! DO j=1,N_SEG+1
	! CUR = (N_ROD-1)*(N_SEG+1) + j
	! CURM= (N_ROD-2)*(N_SEG+1) +j

	! TEMP_P = -1.0 + tanh( (COORDS(CUR,2)-COORDS(CURM,2))/WY )
	! TEMP_P = TEMP_P - (2-COORDS(CUR,2)+COORDS(CURM,2))*( 1.0/WY * 1.0/(COSH((COORDS(CUR,2)-COORDS(CURM,2))/WY))**2 )
	! TEMP_P = TEMP_P*KY*0.5

  !  	GRAD_SPAT(CUR,2)=TEMP_P
  ! ENDDO

  
  !6.2.1) Now convert these to dE/dtheta_i
  DO i=1,N_ROD
	DO j=1,N_SEG-1
		! 1D location of current angle
		CUR=(i-1)*(N_SEG-1)+j
		DO k=j+1,N_SEG
			! 1D location of a node which the current angle affects
			NODE_K=(i-1)*(N_SEG+1)+k+1
			GRAD_A(CUR) = GRAD_A(CUR) - LENGTHS(CUR)*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LENGTHS(CUR)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
		ENDDO		
	ENDDO
  ENDDO


  !6.2.2) Now convert these to dE/dlength_i
    DO i=1,N_ROD
	DO j=1,N_SEG-1
  	! 1D location of current length
		CUR=(i-1)*(N_SEG-1)+j
		DO k=j+1,N_SEG
			! 1D location of a node which the current length affects
			NODE_K=(i-1)*(N_SEG+1)+k+1
			GRAD_L(CUR) = GRAD_L(CUR) + COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
		ENDDO		
	ENDDO
  ENDDO
  !PRINT *, GRAD_L

!   !6.2) Compute the angular and length gradients 
!   !6.2.1) Interior rods
!   !Angles:
!   DO i=2,N_ROD-1
! 	DO j=1,N_SEG-1
! 		! 1D location of current angle
! 		CUR=(i-1)*(N_SEG-1)+j
	
! 		DO k=j+1,N_SEG+1
! 			! 1D location of a node which the current angle affects
! 			NODE_K=(i-1)*(N_SEG+1)+k+1
			
! 			!GRAD(CUR) = GRAD(CUR) + LENGTHS(j)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
!       GRAD_A(CUR) = GRAD_A(CUR) - LENGTHS(CUR)*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LENGTHS(CUR)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
! 		ENDDO		
! 	ENDDO
!   ENDDO
!  !Lengths:
!   DO i=2,N_ROD-1
!     DO j=1,N_SEG-1
!       ! 1D location of current angle
!       CUR=(i-1)*(N_SEG-1)+j
    
!       DO k=j+1,N_SEG+1
!         ! 1D location of a node which the current angle affects
!         NODE_K=(i-1)*(N_SEG+1)+k+1
        
!         !GRAD(CUR) = GRAD(CUR) + LENGTHS(j)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
!         GRAD_L(CUR) = GRAD_L(CUR) + COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
!       ENDDO		
!     ENDDO
!     ENDDO

!   !6.2.2) Exterior rods
!   !Angles:
! 	DO j=1,N_SEG-1
! 		! 1D location of current angle
! 		CUR=j
	
! 		DO k=j+1,N_SEG+1
! 			!First rod
! 			NODE_K=k+1
! 			GRAD_A(CUR) = GRAD_A(CUR) - LENGTHS(CUR)*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1)+ LENGTHS(CUR)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
! 		ENDDO
! 	ENDDO
  
!   DO j=1,N_SEG-1
! 		! 1D location of current angle
! 		CUR=(N_ROD-1)*(N_SEG-1)+j
	
! 		DO k=j+1,N_SEG+1
! 			!First rod
! 			NODE_K=(N_ROD-1)*(N_SEG+1)+k+1
! 			GRAD_A(CUR) = GRAD_A(CUR) - LENGTHS(CUR)*SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + LENGTHS(CUR)*COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
! 		ENDDO
! 	ENDDO
!   !Lengths:
!   DO j=1,N_SEG-1
! 		! 1D location of current length
! 		CUR=j
	
! 		DO k=j+1,N_SEG+1
! 			!First rod
! 			NODE_K=k+1
! 			GRAD_L(CUR) = GRAD_L(CUR) + COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,1)+ SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
! 		ENDDO
! 	ENDDO
  
!   DO j=1,N_SEG-1
! 		! 1D location of current length
! 		CUR=(N_ROD-1)*(N_SEG-1)+j
	
! 		DO k=j+1,N_SEG+1
! 			!First rod
! 			NODE_K=(N_ROD-1)*(N_SEG+1)+k+1
! 			GRAD_A(CUR) = GRAD_A(CUR) +COS(ANGLES(CUR))*GRAD_SPAT(NODE_K,1) + SIN(ANGLES(CUR))*GRAD_SPAT(NODE_K,2)
! 		ENDDO
! 	ENDDO

! !   !7) Energy penalising an angle geeing above pi
  
! !   !DO i=1,N_ROD
! ! 	!DO j=1,N_SEG-1
! ! 	!	CUR=(i-1)*(N_SEG-1)+j
! ! 	!	ETH = ETH + KY*TANH((ABS(ANGLES(CUR))-3.14159265359)/0.3);
! !   ! GRAD(CUR) = GRAD(CUR)+
! ! 	!ENDDO
! !   !ENDDO

!   !8) Energy penalising length being 0 -- very soft constraint... :(
!HARMONIC-----------
  DO i=1,N_ROD
    DO j=1,N_SEG+1
      CUR = (i-1)*(N_SEG-1) + j
      GRAD_L(CUR) = GRAD_L(CUR) + KL/2*(2*(LENGTHS(CUR)-LS_0))
    ENDDO
 ENDDO
! !LJ------------------
! DO i=1,N_ROD
!   DO j=1, N_SEG+1
!     CUR = (i-1)*(N_SEG-1) + j
!     GRAD_L(CUR)=GRAD_L(CUR)+(KL**2)*((-12*LS_0**12)/(LENGTHS(CUR)**13)+(6*LS_0**6)/(LENGTHS(CUR)**7))

!   ENDDO
! ENDDO


!  ---------------------------
!open(2, file='grad.out')
!write(2,'(F20.10)') GRAD
!close(2)
!stop
END SUBROUTINE COMPUTE_GRAD
!------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------
!Test the gradient vs a finite difference approximation
SUBROUTINE TEST_GRADIENT(G,X)

! IMPLICIT NONE
! INTEGER :: I
! DOUBLE PRECISION :: DT, A_0, A_P, A_M, E_P, E_M, TEST_GRAD(N), GRAD_DIFF(N), G(N)

! DT=0.000001

! DO I=1,N
! 	A_0 = ANGLES(I)
! 	A_P = ANGLES(I) + DT
! 	A_M = ANGLES(I) - DT

! 	E=0.0

! 	ANGLES(I) = A_P
! 	CALL UNPACK_COORDINATES()
! 	CALL COMPUTE_ENERGY(E)
! 	E_P=E
! 	ANGLES(I)=A_M
! 	CALL UNPACK_COORDINATES()
! 	CALL COMPUTE_ENERGY(E)
! 	E_M=E

! 	TEST_GRAD(I) = (E_P-E_M)/(2*DT)*CM1	
	
! 	ANGLES(I) = A_0

! 	!PRINT *, G(I), TEST_GRAD(I), G(I)-TEST_GRAD(I)
!  ! IF (ABS((G(I)-TEST_GRAD(I))/G(I))>1.0D-5) THEN
! 	!PRINT *, 'TEST GRAD FAILED AT I=', I
! 	!ENDIF
! ENDDO

! GRAD_DIFF = (G-TEST_GRAD)/G



! OPEN(1,FILE='grad.diff')
! WRITE(1,'(F20.10)') GRAD_DIFF
! CLOSE(1)

! OPEN(1,FILE='grad.orig')
! WRITE(1,'(F20.10)') G
! CLOSE(1)

! OPEN(1,FILE='grad.test')
! WRITE(1,'(F20.10)') TEST_GRAD
! CLOSE(1)

!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: I
DOUBLE PRECISION :: DT, A_0, A_P, A_M, L_0, L_P, L_M, E_P, E_M, TEST_GRAD(N), GRAD_DIFF(N), G(N), X(N)
!PRINT *, "X", X
ANGLES(:)=X(1:N/2)
LENGTHS(:)=X(N/2+1:N)
!PRINT *, "LENGTHS(1)", LENGTHS(1)
!PRINT *, "X", X
DT=0.000001
 CALL COMPUTE_GRAD()
 G(1:N/2)=GRAD_A
 G(N/2+1:N)=GRAD_L
 !PRINT *, GRAD_L
 
DO I=1,N
 !PRINT *, "X diff", X(I)

 !PRINT *, "X diff", X(1)
  !PRINT *, "diff, ANGLES= " , ANGLES
IF (I<(N/2)+1) THEN
!PRINT *, "X(",I,") = ", X(I)
	A_0 = ANGLES(I)
	A_P = ANGLES(I) + DT
	A_M = ANGLES(I) - DT
  E=0.0
  ANGLES(I) = A_P
  !PRINT *, "ANGLES I", ANGLES
	CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
  E_P=E
  !PRINT *, E_P
  ANGLES(I)=A_M
  CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
  E_M=E
  !PRINT *, E_M
  TEST_GRAD(I) = (E_P-E_M)/(2*DT)
  !Print *, E_P
  !Print *, E_M

  GRAD_DIFF(I) = (G(I)-TEST_GRAD(I))/G(I)
  !PRINT *, "GRAD TEST(", I, ")", TEST_GRAD(I)
  !PRINT *, "IN", TEST_GRAD
  ANGLES(I) = A_0
ENDIF

IF (I>(N/2)) THEN
 !PRINT *, "LENGTHS(", I-N/2, ")", LENGTHS(I-N/2)
  L_0 = LENGTHS(I-N/2)
	L_P = LENGTHS(I-N/2) + DT
	L_M = LENGTHS(I-N/2) - DT
  !PRINT *, "STARTING"
  E=0.0
  LENGTHS(I-N/2) = L_P
	CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
  !PRINT*, "CALLED"
  !PRINT *, E
  E_P=E
  !PRINT *, "EP", E_P
  LENGTHS(I-N/2)=L_M
  !PRINT *, "LENGTHS(",I-N/2,")", LENGTHS(I-N/2)
  CALL UNPACK_COORDINATES()
	CALL COMPUTE_ENERGY(E)
  E_M=E
    !PRINT *, "EM", E_M
  TEST_GRAD(I) = (E_P-E_M)/(2*DT)
   !  PRINT *, "G", G(I)
   
  GRAD_DIFF(I) = (G(I)-TEST_GRAD(I))/G(I)
   !PRINT *, "gdiff", GRAD_DIFF(I)
   !PRINT *, "gtest", TEST_GRAD(I)
  LENGTHS(I-N/2) = L_0

  ! PRINT *, "GRAD TEST(", I, ")", TEST_GRAD(I)
   !PRINT *, "E_P", E_P
     !PRINT *, "E_M", E_M
   ! PRINT *, "L_P", L_P
    ! PRINT *, "L_M", L_M
    
ENDIF

	!PRINT *, G(I), TEST_GRAD(I), G(I)-TEST_GRAD(I)
	!IF (ABS(G(I)-TEST_GRAD(I))>1.0D-5) THEN
	!	PRINT *, 'TEST GRAD FAILED AT I=', I
	!ENDIF
   !PRINT *, "grad diff computed"
   !PRINT *, TEST_GRAD

ENDDO





OPEN(1,FILE='grad.diff')
WRITE(1,'(F20.10)') GRAD_DIFF
CLOSE(1)
!PRINT *, "G", G
OPEN(1,FILE='grad.orig')
WRITE(1,'(F25.10)') G
CLOSE(1)
!PRINT*, "OUT", TEST_GRAD
OPEN(1,FILE='grad.test')
WRITE(1,'(F30.10)') TEST_GRAD
CLOSE(1)

STOP

END SUBROUTINE TEST_GRADIENT


END MODULE potential
