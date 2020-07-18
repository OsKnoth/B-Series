MODULE MISBI_Mod
  USE BiColourTree_Mod
  USE BSeries_Mod
  IMPLICIT NONE

  TYPE Koeff_T 
    REAL(8), ALLOCATABLE :: a(:)
  END TYPE Koeff_T 
  TYPE MISBIMethod_T
    INTEGER :: nStage
    INTEGER :: nPhi
    INTEGER, ALLOCATABLE :: sPhi(:)
    TYPE(Koeff_T), ALLOCATABLE :: a(:,:)
    REAL(8), ALLOCATABLE :: d(:,:)
    REAL(8), ALLOCATABLE :: g(:,:)
    REAL(8), ALLOCATABLE :: c(:)
    REAL(8), ALLOCATABLE :: dt(:)
    REAL(8), ALLOCATABLE :: cw(:)
  END TYPE MISBIMethod_T
  LOGICAL :: dInsert=.TRUE.
  LOGICAL :: gInsert=.TRUE.
  INTEGER :: nStageSym=4
  INTEGER :: nPhiSym=4


CONTAINS  

SUBROUTINE MISBIBSeriesCompose1(PhiMISBI,MISBI,pMax)
  TYPE(BSeries_T), ALLOCATABLE :: PhiMISBI(:)
  TYPE(MISBIMethod_T) :: MISBI
  INTEGER :: pMax

  INTEGER :: i,is,j,js,k,ns,p
  INTEGER :: NumberLoc,OrderLoc
  TYPE(BSeries_T), ALLOCATABLE :: Eta(:)
  TYPE(BPSeries_T), ALLOCATABLE :: Z(:)
  TYPE(BSeries_T), ALLOCATABLE :: PhiMISBIb(:)
  TYPE(BSeries_T) :: TempB
  TYPE(BPSeries_T) :: TempW
  TYPE(Tree_T), POINTER :: Tree
  TYPE(TreeP_T), POINTER :: TreeP
  TYPE(Tree_T), POINTER :: Child
  TYPE(PKoeff_T) :: PowerOne

  !CALL PowerOnePolynom(PowerOne)
  nS=MISBI%nStage
  ALLOCATE(PhiMISBI(nS+1))
  ALLOCATE(Eta(nS+1))
  ALLOCATE(Z(nS+1))
  ALLOCATE(PhiMISBIb(nS+1))
  DO iS=1,nS+1
    CALL AllocateBSeries(Eta(iS),pMax)
    CALL AllocateBPSeries(Z(iS),pMax)
    CALL AllocateBSeries(PhiMISBI(iS),pMax)
    Eta(iS)%Phi(0)%a(1)=1.0d0
    Z(iS)%Phi(0)%a(1)%Koeff=0.0d0
    Z(iS)%Phi(0)%a(1)%Koeff(0)=1.0d0
    Z(iS)%Phi(0)%a(1)%Grad=0
  END DO  
  CALL AllocateBSeries(TempB,pMax)
  CALL AllocateBPSeries(TempW,pMax)
  DO p=1,pMax
    DO iS=1,ns+1
      DO jS=1,ns
        TreeP=>ListTree(p)%T
        DO i=1,ListTree(p)%LenListTree
          Tree=>TreeP%TP
          CALL UnitPolynom(TempW%Phi(p)%a(i))
          DO j=1,Tree%NumberChilds
            Child=>Tree%Childs(j)%TP
            OrderLoc=Child%Order
            NumberLoc=Child%NumberOrder
            DO k=1,Tree%OrderChilds(j)
              IF (Child%Colour=='w') THEN
                CALL MultPolynom(TempW%Phi(p)%a(i),TempW%Phi(p)%a(i),Z(iS)%Phi(OrderLoc)%a(NumberLoc))
              ELSE
                CALL MultScalarPolynom(TemPW%Phi(p)%a(i),Eta(jS)%Phi(OrderLoc)%a(NumberLoc))
              END IF  
            END DO  
          END DO  
          DO k=1,MISBI%sPhi(iS)
            CALL AddPolynomScalar(Z(iS)%Phi(p)%a(i),Z(iS)%Phi(p)%a(i),MISBI%A(is,js)%a(k),TempW%Phi(p)%a(i))
            IF (k<MISBI%nPhi) THEN
              CALL MultPolynom(TempW%Phi(p)%a(i),TempW%Phi(p)%a(i),PowerOne)
            END IF
          END DO
          TreeP=>TreeP%Next
        END DO
      END DO  
    END DO  
    DO iS=1,ns+1
      TreeP=>ListTree(p)%T
      DO i=1,ListTree(p)%LenListTree
        Tree=>TreeP%TP
        DO jS=1,nS
          Z(iS)%Phi(p)%a(i)%Koeff(0)=Z(iS)%Phi(p)%a(i)%Koeff(0)+MISBI%g(is,js)*(Eta(jS)%Phi(p)%a(i)-Eta(1)%Phi(p)%a(i))
        END DO
        CALL IntegratePolynom(Z(iS)%Phi(p)%a(i),Z(iS)%Phi(p)%a(i))
        Z(iS)%Phi(p)%a(i)%Koeff(0)=Eta(1)%Phi(p)%a(i)
        DO jS=1,ns
          Z(iS)%Phi(p)%a(i)%Koeff(0)=Z(iS)%Phi(p)%a(i)%Koeff(0)+MISBI%d(is,js)*(Eta(jS)%Phi(p)%a(i)-Eta(1)%Phi(p)%a(i))
        END DO
        Eta(iS)%Phi(p)%a(i)=ValuePolynom(Z(iS)%Phi(p)%a(i),1.0d0)
        TreeP=>TreeP%Next
      END DO
    END DO  
  END DO  
  DO p=0,pMax
    DO iS=1,nS+1
      PhiMISBI(iS)%Phi(p)%a(:)=Eta(iS)%Phi(p)%a(:)
    END DO
  END DO  

END SUBROUTINE MISBIBSeriesCompose1

SUBROUTINE MISBIBSeriesCompose(PhiMISBI,PhiPMISBI,MISBI,pMax)
  TYPE(BSeries_T), ALLOCATABLE :: PhiMISBI(:)
  TYPE(BPSeries_T), ALLOCATABLE :: PhiPMISBI(:)
  TYPE(MISBIMethod_T) :: MISBI
  INTEGER :: pMax

  INTEGER :: i,is,j,js,k,ns,p
  INTEGER :: NumberLoc,OrderLoc
  TYPE(BSeries_T), ALLOCATABLE :: Eta(:)
  TYPE(BPSeries_T), ALLOCATABLE :: Z(:)
  TYPE(BSeries_T), ALLOCATABLE :: PhiMISBIb(:)
  TYPE(BSeries_T) :: TempB
  TYPE(BPSeries_T) :: TempW
  TYPE(Tree_T), POINTER :: Tree
  TYPE(TreeP_T), POINTER :: TreeP
  TYPE(Tree_T), POINTER :: Child
  TYPE(PKoeff_T) :: PowerOne(1:MISBI%nPhi)
  CHARACTER :: isW,jsW,kW,pW
  CHARACTER(30) :: String

  INTEGER :: l

  CALL PowerOnePolynom(PowerOne)
  nS=MISBI%nStage
  ALLOCATE(PhiMISBI(nS+1))
  ALLOCATE(PhiPMISBI(nS+1))
  ALLOCATE(Eta(nS+1))
  ALLOCATE(Z(nS+1))
  ALLOCATE(PhiMISBIb(nS+1))
  DO iS=1,nS+1
    CALL AllocateBSeries(Eta(iS),pMax)
    CALL AllocateBPSeries(Z(iS),pMax)
    CALL AllocateBSeries(PhiMISBI(iS),pMax)
    CALL AllocateBPSeries(PhiPMISBI(iS),pMax)
    Eta(iS)%Phi(0)%a(1)=1.0d0
    Z(iS)%Phi(0)%a(1)%Koeff=0.0d0
    Z(iS)%Phi(0)%a(1)%Koeff(0)=1.0d0
    Z(iS)%Phi(0)%a(1)%Grad=0
  END DO  
  CALL AllocateBSeries(TempB,pMax)
  CALL AllocateBPSeries(TempW,pMax)
  DO p=1,pMax
    DO iS=1,ns+1
      DO jS=1,ns
        TreeP=>ListTree(p)%T
        DO i=1,ListTree(p)%LenListTree
          Tree=>TreeP%TP
          CALL UnitPolynom(TempW%Phi(p)%a(i))
          DO j=1,Tree%NumberChilds
            Child=>Tree%Childs(j)%TP
            OrderLoc=Child%Order
            NumberLoc=Child%NumberOrder
            DO k=1,Tree%OrderChilds(j)
              IF (Child%Colour=='w') THEN
                CALL MultPolynom(TempW%Phi(p)%a(i),TempW%Phi(p)%a(i),Z(iS)%Phi(OrderLoc)%a(NumberLoc))
              ELSE  
                CALL MultScalarPolynom(TempW%Phi(p)%a(i),Eta(jS)%Phi(OrderLoc)%a(NumberLoc))
              END IF  
            END DO  
          END DO  
          DO k=1,MISBI%sPhi(iS)
            CALL AddPolynomScalar(Z(iS)%Phi(p)%a(i),Z(iS)%Phi(p)%a(i),MISBI%A(is,js)%a(k),TempW%Phi(p)%a(i))
            IF (k<MISBI%sPhi(iS)) THEN
              CALL MultPolynom(TempW%Phi(p)%a(i),TempW%Phi(p)%a(i),PowerOne(k))
            END IF
          END DO
          TreeP=>TreeP%Next
        END DO
      END DO  
    END DO  
    DO iS=1,ns+1
      TreeP=>ListTree(p)%T
      DO i=1,ListTree(p)%LenListTree
        Tree=>TreeP%TP
        DO jS=1,ns
          Z(iS)%Phi(p)%a(i)%Koeff(0)=Z(iS)%Phi(p)%a(i)%Koeff(0)+MISBI%g(is,js)*(Eta(jS)%Phi(p)%a(i)-Eta(1)%Phi(p)%a(i))
        END DO
        CALL IntegratePolynom(Z(iS)%Phi(p)%a(i),Z(iS)%Phi(p)%a(i))
        Z(iS)%Phi(p)%a(i)%Koeff(0)=0.0d0
        DO jS=1,ns
          Z(iS)%Phi(p)%a(i)%Koeff(0)=Z(iS)%Phi(p)%a(i)%Koeff(0)+MISBI%d(is,js)*(Eta(jS)%Phi(p)%a(i)-Eta(1)%Phi(p)%a(i))
        END DO
        Eta(iS)%Phi(p)%a(i)=ValuePolynom(Z(iS)%Phi(p)%a(i),1.0d0)
        TreeP=>TreeP%Next
      END DO
    END DO  
  END DO  
  DO p=0,pMax
    DO iS=1,nS+1
      PhiMISBI(iS)%Phi(p)%a(:)=Eta(iS)%Phi(p)%a(:)
      DO k=LBOUND(PhiMISBI(iS)%Phi(p)%a,1),UBOUND(PhiMISBI(iS)%Phi(p)%a,1)
        PhiPMISBI(iS)%Phi(p)%a(k)=Z(iS)%Phi(p)%a(k)
      END DO  
    END DO
  END DO  

END SUBROUTINE MISBIBSeriesCompose
SUBROUTINE MISBIMethod(Method,MISBI,pMax)
  CHARACTER(*) :: Method
  TYPE(MISBIMethod_T) :: MISBI
  INTEGER :: pMax

  INTEGER :: i,j,k
  REAL(8) :: Fac
  REAL(8) :: delta
  REAL(8) :: Rand
  REAL(8) :: SumRow
  REAL(8), ALLOCATABLE :: A(:)
  REAL(8), ALLOCATABLE :: B(:)

  WRITe(*,*) 'Method ',Method
  SELECT CASE(Method)
    CASE ('ExSym')
      MISBI%nStage=nStageSym
      MISBI%nPhi=nPhiSym
      CALL AllocateMISBI
      MISBI%sPhi=nPhiSym
      CALL RANDOM_SEED()
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          CALL RANDOM_NUMBER(Rand)
          IF (gInsert) THEN
            MISBI%g(i,j)=Rand
            MISBI%g(i,j)=0.0d0
          END IF  
          IF (dInsert) THEN
            MISBI%d(i,j)=Rand
            MISBI%d(i,j)=0.0d0
          END IF  
          DO k=1,nPhiSym
            CALL RANDOM_NUMBER(Rand)
            MISBI%a(i,j)%a(k)=Rand
            MISBI%a(i,j)%a(k)=5.0d0
          END DO
        END DO
      END DO
      DO i=2,MISBI%nStage+1
        MISBI%c(i)=0.0d0
        DO j=1,i-1
          MISBI%c(i)=MISBI%c(i)+MISBI%a(i,j)%a(1)
        END DO  
        DO k=2,nPhiSym
          SumRow=0.0d0
          DO j=1,i-1
            SumRow=SumRow+MISBI%a(i,j)%a(k)
          END DO  
          SumRow=SumRow/(i-1)
          DO j=1,i-1
            MISBI%a(i,j)%a(k)=MISBI%a(i,j)%a(k)-SumRow
          END DO  
        END DO
      END DO
      MISBI%c(MISBI%nStage+1)=1.0d0
      SumRow=0.0d0
      DO i=1,MISBI%nStage
        SumRow=SumRow+MISBI%a(nStageSym+1,i)%a(1)
      END DO  
      DO i=1,MISBI%nStage
        MISBI%a(MISBI%nStage+1,i)%a(1)=MISBI%a(MISBI%nStage+1,i)%a(1)/SumRow
      END DO  
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE('ExEuler') 
      MISBI%nStage=1
      MISBI%nPhi=1
      MISBI%sPhi=MISBI%nPhi
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=1.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0
    CASE('ExRK2a') 
      MISBI%nStage=2
      MISBI%nPhi=1
      MISBI%sPhi=MISBI%nPhi
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=0.5d0
      MISBI%a(3,1)%a(1)=-0.5d0
      MISBI%a(3,2)%a(1)=1.0d0

      MISBI%d(3,2)=1.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=0.5d0
      MISBI%dt(3)=0.5d0
    CASE('ExRK2a1') 
      MISBI%nStage=2
      MISBI%nPhi=1
      MISBI%sPhi=MISBI%nPhi
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=0.5d0
      MISBI%a(3,2)%a(1)=1.0d0

      MISBI%d(3,2)=0.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=0.5d0
      MISBI%dt(3)=1.0d0
    CASE('ExRK2b') 
      MISBI%nStage=2
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=1.0d0
      MISBI%a(3,1)%a(1)=0.5d0
      MISBI%a(3,2)%a(1)=0.5d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0
      MISBI%dt(3)=1.0d0
    CASE('ExHeun')
      MISBI%nStage=3
      MISBI%nPhi=3
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=1.0d0/3.0d0
      MISBI%sPhi(2)=1
      MISBI%a(3,2)%a(1)=2.0d0/3.0d0
      MISBI%sPhi(3)=2
      MISBI%a(4,1)%a(1)=1.0d0/4.0d0
      MISBI%a(4,3)%a(1)=3.0d0/4.0d0
      MISBI%sPhi(4)=3
      MISBI%c(1)=0.0d0
      MISBI%c(2)=1.0d0/3.0d0
      MISBI%c(3)=2.0d0/3.0d0
      MISBI%c(4)=1.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0
      MISBI%dt(3)=1.0d0
    CASE('ExRK3Nair')  
      MISBI%nStage=3
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=1.0d0/3.0d0
      MISBI%a(3,2)%a(1)=1.0d0/6.0d0
      MISBI%a(4,3)%a(1)=1.0d0
      MISBI%d(3,2)=1.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0/3.0d0
      MISBI%dt(3)=1.0d0/6.0d0
      MISBI%dt(4)=1.0d0
    CASE('ExRK3')  
      MISBI%nStage=3
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(2,1)%a(1)=1.0d0/3.0d0
      MISBI%a(3,2)%a(1)=1.0d0/2.0d0
      MISBI%a(4,3)%a(1)=1.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0/3.0d0
      MISBI%dt(3)=1.0d0/2.0d0
      MISBI%dt(4)=1.0d0
    CASE('ExRK4')
      MISBI%nStage=4
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=0.5d0
      MISBI%a(3,2)%a(1)=0.5d0
      MISBI%a(4,3)%a(1)=1.0d0
      MISBI%a(5,1)%a(1)=1.0d0/6.0d0
      MISBI%a(5,2)%a(1)=2.0d0/6.0d0
      MISBI%a(5,3)%a(1)=2.0d0/6.0d0
      MISBI%a(5,4)%a(1)=1.0d0/6.0d0
      MISBI%c(1)=0.0d0
      MISBI%c(2)=0.5d0
      MISBI%c(3)=0.5d0
      MISBI%c(4)=1.0d0
      MISBI%c(5)=1.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=0.5d0
      MISBI%dt(3)=0.5d0
      MISBI%dt(4)=1.0d0
      MISBI%dt(5)=1.0d0
    CASE('ExRK2a2') 
      MISBI%nStage=2
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%a(2,1)%a(1)=0.0d0
      MISBI%a(3,2)%a(1)=0.0d0

      MISBI%d(3,2)=0.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0
      MISBI%dt(3)=1.0d0
    CASE('Implicit2') 
      MISBI%nStage=2
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(1,1)%a(1)=5.0d0/12.0d0
      MISBI%a(1,2)%a(1)=-1.0d0/12.0d0
      MISBI%a(2,1)%a(1)=3.0d0/4.0d0
      MISBI%a(2,2)%a(1)=1.0d0/4.0d0
      MISBI%a(3,1)%a(1)=3.0d0/4.0d0
      MISBI%a(3,2)%a(1)=1.0d0/4.0d0


      MISBI%dt(1)=1.0d0/3.0d0
      MISBI%dt(2)=1.0d0
      MISBI%dt(3)=1.0d0
    CASE('Gauss2')  
      MISBI%nStage=2
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(1,1)%a(1)=1.0d0/4.0d0
      MISBI%a(1,2)%a(1)=(3.0d0-2.0d0*SQRT(3.0d0))/12.0d0
      MISBI%a(2,1)%a(1)=(3.0d0+2.0d0*SQRT(3.0d0))/12.0d0
      MISBI%a(2,2)%a(1)=1.0d0/4.0d0
      MISBI%a(3,1)%a(1)=1.0d0/2.0d0
      MISBI%a(3,2)%a(1)=1.0d0/2.0d0
      MISBI%dt(1)=(3-SQRT(3.0d0))/6.0d0
      MISBI%dt(2)=(3+SQRT(3.0d0))/6.0d0
      MISBI%dt(3)=1.0d0
    CASE ('KWRK43')
      WRITE(*,*) 'Methode KWRK43'
      MISBI%nStage=4
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(2,1)%a(1)=1.0d0/2.0d0
      MISBI%a(3,1)%a(1)=-1.0d0/6.0d0
      MISBI%a(3,2)%a(1)=2.0d0/3.0d0
      MISBI%a(4,1)%a(1)=1.0d0/3.0d0
      MISBI%a(4,2)%a(1)=-1.0d0/3.0d0
      MISBI%a(4,3)%a(1)=1.0d0
      MISBI%a(5,1)%a(1)=1.0d0/6.0d0
      MISBI%a(5,2)%a(1)=1.0d0/3.0d0
      MISBI%a(5,3)%a(1)=1.0d0/3.0d0
      MISBI%a(5,4)%a(1)=1.0d0/6.0d0
      DO j=1,MISBI%nStage
        MISBI%a(5,j)%a=MISBI%a(5,j)%a-MISBI%a(4,j)%a
        MISBI%a(4,j)%a=MISBI%a(4,j)%a-MISBI%a(3,j)%a
        MISBI%a(3,j)%a=MISBI%a(3,j)%a-MISBI%a(2,j)%a
      END DO
      MISBI%c=0.0d0
      MISBI%c(2)=1.0d0/2.0d0
      MISBI%c(3)=1.0d0/2.0d0
      MISBI%c(4)=1.0d0
      MISBI%c(5)=1.0d0
      MISBI%d=0.0d0
      MISBI%d(2,1)=1.0d0
      MISBI%d(3,2)=1.0d0
      MISBI%d(4,3)=1.0d0
      MISBI%d(5,4)=1.0d0
      MISBI%g=0.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=0.5d0
      MISBI%dt(3)=0.0d0
      MISBI%dt(4)=0.5d0
      MISBI%dt(5)=0.0d0
    CASE ('MISBI4_4')
      MISBI%nStage=4
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(2,1)%a(1)=  0.38758444641450318     
      MISBI%a(3,1)%a(1)=  -2.5318448354142823E-002
      MISBI%a(3,2)%a(1)=  0.38668943087310403     
      MISBI%a(4,1)%a(1)=  0.20899983523553325     
      MISBI%a(4,2)%a(1)= -0.45856648476371231     
      MISBI%a(4,3)%a(1)=  0.43423187573425748     
      MISBI%a(5,1)%a(1)= -0.10048822195663100     
      MISBI%a(5,2)%a(1)= -0.46186171956333327     
      MISBI%a(5,3)%a(1)=  0.83045062122462809     
      MISBI%a(5,4)%a(1)=  0.27014914900250392     
      MISBI%c(2)=  0.38758444641450318     
      MISBI%c(3)=  0.61521685655017821     
      MISBI%c(4)=  0.23254717315441453     
      MISBI%c(5)=   1.0000000000000002     
      MISBI%d(3,2)=  0.52349249922385610     
      MISBI%d(4,2)=   1.1683374366893629     
      MISBI%d(4,3)= -0.75762080241712637     
      MISBI%d(5,2)=  -3.6477233846797109E-002
      MISBI%d(5,3)=  0.56936148730740477     
      MISBI%d(5,4)=  0.47746263002599681     
      MISBI%g(3,2)=  0.13145089796226542     
      MISBI%g(4,2)= -0.36855857648747881     
      MISBI%g(4,3)=  0.33159232636600550     
      MISBI%g(5,2)=  -6.5767130537473045E-002
      MISBI%g(5,3)=   4.0591093109036858E-002
      MISBI%g(5,4)=   6.4902111640806712E-002
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('MISBI54_IIB1')
      MISBI%nStage=5
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1

      MISBI%a(2,1)%a(1)=  0.219579314792533d0
      MISBI%a(3,1)%a(1)= -0.032864918414060d0
      MISBI%a(3,2)%a(1)=  0.634699918767414d0
      MISBI%a(4,1)%a(1)= -0.241761887431829d0
      MISBI%a(4,2)%a(1)= -0.120631540663984d0
      MISBI%a(4,3)%a(1)=  0.374686620841487d0
      MISBI%a(5,1)%a(1)= -0.058474324094343d0
      MISBI%a(5,2)%a(1)=  0.351217252190521d0
      MISBI%a(5,3)%a(1)=  0.309657030167295d0
      MISBI%a(5,4)%a(1)=  0.168604799122988d0
      MISBI%a(6,1)%a(1)= -0.056205055946158d0
      MISBI%a(6,2)%a(1)= -0.068390330952311d0
      MISBI%a(6,3)%a(1)= -0.086209210260269d0
      MISBI%a(6,4)%a(1)=  0.034904705602768d0
      MISBI%a(6,5)%a(1)=  0.448964988009822d0
          

      MISBI%c(2)=  0.219579314792533d0
      MISBI%c(3)=  0.618448183300698d0
      MISBI%c(4)=  0.530092519575536d0
      MISBI%c(5)=  0.808136457883832d0
      MISBI%c(6)=  1.000000000000000d0
      
!     MISBI%d(2,1)= -0.056843003311023d0
      MISBI%d(2,1)=  0.0d0
      MISBI%d(3,1)=  0.071035715986068d0
      MISBI%d(3,1)=  0.0d0
      MISBI%d(3,2)=  0.050143439731979d0
      MISBI%d(4,1)=  0.021491523917140d0
      MISBI%d(4,1)=  0.0d0
      MISBI%d(4,2)=  0.287530720188756d0
      MISBI%d(4,3)=  0.239030810792355d0
      MISBI%d(5,1)=  0.027558616966568d0
      MISBI%d(5,1)=  0.0d0
      MISBI%d(5,2)=  0.382675659910308d0
      MISBI%d(5,3)=  0.177185696263246d0
      MISBI%d(5,4)= -0.314894383613333d0
      MISBI%d(6,1)=  0.065158401284120d0
      MISBI%d(6,1)=  0.0d0
      MISBI%d(6,2)=  0.079591607322196d0
      MISBI%d(6,3)=  0.459806401597571d0
      MISBI%d(6,4)=  0.086725275506356d0
      MISBI%d(6,5)=  0.439945196292364d0

      MISBI%g(2,1)=  0.168489083931286d0
      MISBI%g(2,1)=  0.0d0
      MISBI%g(3,1)= -0.025097850341834d0
      MISBI%g(3,1)=  0.0d0
      MISBI%g(3,2)=  0.025515704040468d0
      MISBI%g(4,1)=  0.106139356407192d0
      MISBI%g(4,1)=  0.0d0
      MISBI%g(4,2)=  0.264445452990869d0
      MISBI%g(4,3)=  0.402246482358727d0
      MISBI%g(5,1)= -0.031464053194458d0
      MISBI%g(5,1)=  0.0d0
      MISBI%g(5,2)= -0.068258296801680d0
      MISBI%g(5,3)=  0.027558616966568d0
      MISBI%g(5,4)=  0.015830368641068d0
      MISBI%g(6,1)=  0.150547662349659d0
      MISBI%g(6,1)=  0.0d0
      MISBI%g(6,2)=  0.088610905686011d0
      MISBI%g(6,3)=  0.067880982803316d0
      MISBI%g(6,4)= -0.297416190393485d0
      MISBI%g(6,5)=  0.148246909195494d0
    
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('MRI-GARK-ERK33')
    ! Sandu
      WRITE(*,*) 'Method: MRI-GARK-ERK33'
      MISBI%nStage=3
      MISBI%nPhi=2
      CALL AllocateMISBI
    ! MISBI%sPhi(2)=1
    ! MISBI%sPhi(3)=2
    ! MISBI%sPhi(4)=2
      delta=-0.5d0
      MISBI%a(2,1)%a(1)=1.0d0/3.0d0
      MISBI%a(3,1)%a(1)=(-6.0d0*delta-7.0d0)/12.0d0
      MISBI%a(3,2)%a(1)=(6.0d0*delta+11.0d0)/12.0d0
      MISBI%a(4,2)%a(1)=(6.0d0*delta-5.0d0)/12.0d0
      MISBI%a(4,3)%a(1)=(3.0d0-2.0d0*delta)/4.0d0

      MISBI%a(3,1)%a(2)=(2.0d0*delta+1.0d0)/2.0d0
      MISBI%a(3,2)%a(2)=-(2.0d0*delta+1.0d0)/2.0d0
      MISBI%a(4,1)%a(2)=0.5d0
      MISBI%a(4,2)%a(2)=-(2.0d0*delta+1.0d0)/2.0d0
      MISBI%a(4,3)%a(2)=delta

      MISBI%c=0.0d0
      MISBI%c(2)=1.0d0/3.0d0
      MISBI%c(3)=2.0d0/2.0d0
      MISBI%c(4)=1.0d0
      MISBI%d=0.0d0
      MISBI%d(2,1)=1.0d0
      MISBI%d(3,2)=1.0d0
      MISBI%d(4,3)=1.0d0
      MISBI%g=0.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0/3.0d0
      MISBI%dt(3)=1.0d0/3.0d0
      MISBI%dt(4)=1.0d0/3.0d0
    CASE ('BWRRK33')
      !J. Williamson, Low-storage runge-kutta schemes, Journal of
      !Computational Physics 35 (1) (1980) 48 { 56. doi:https:
      !//doi.org/10.1016/0021-9991(80)90033-9.
      !URL http://www.sciencedirect.com/science/article/pii/0021999180900339
      WRITE(*,*) 'Method BWRRK33 '
      MISBI%nStage=3
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(1)=0.0d0
      A(2)=-0.637694471842202d0
      A(3)=-1.306647717737108d0
      B(1)=0.457379997569388d0
      B(2)=0.925296410920922d0
      B(3)=0.393813594675071d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('CKRK54')
      ! 5-stage fourth-order CKRK54 of
      ! Carpenter, Kennedy, Technical Memorandum NASA-TM-109112 (1994)
      WRITE(*,*) 'Method CKRK54 '
      MISBI%nStage=5
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-567301805773d0/1357537059087d0
      A(3)=-2404267990393d0/2016746695238d0
      A(4)=-3550918686646d0/2091501179385d0
      A(5)=-1275806237668d0/842570457699d0
      B(1)=1432997174477d0/9575080441755d0
      B(2)=5161836677717d0/13612068292357d0
      B(3)=1720146321549d0/2090206949498d0
      B(4)=3134564353537d0/4481467310338d0
      B(5)=2277821191437d0/14882151754819d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO


    CASE ('NDBRKF124')
      ! 12-stage fourth-order NDBRKF124 of
      ! Niegemann, Diehl, Busch, J. Comp. Phys. 231 (2012) 364-372
      WRITE(*,*) 'Method NDBRKF124 '
      MISBI%nStage=12
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.0923311242368072d0
      A(3)=-0.9441056581158819d0
      A(4)=-4.3271273247576394d0
      A(5)=-2.1557771329026072d0
      A(6)=-0.9770727190189062d0
      A(7)=-0.7581835342571139d0
      A(8)=-1.7977525470825499d0
      A(9)=-2.6915667972700770d0
      A(10)=-4.6466798960268143d0
      A(11)=-0.1539613783825189d0
      A(12)=-0.5943293901830616d0
      B(1)=0.0650008435125904d0
      B(2)=0.0161459902249842d0
      B(3)=0.5758627178358159d0
      B(4)=0.1649758848361671d0
      B(5)=0.3934619494248182d0
      B(6)=0.0443509641602719d0
      B(7)=0.2074504268408778d0
      B(8)=0.6914247433015102d0
      B(9)=0.3766646883450449d0
      B(10)=0.0757190350155483d0
      B(11)=0.2027862031054088d0
      B(12)=0.2167029365631842d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE ('NDBRKF134')
      ! 13-stage fourth-order NDBRKF134 of
      ! Niegemann, Diehl, Busch, J. Comp. Phys. 231 (2012) 364-372
      WRITE(*,*) 'Method NDBRKF134 '
      MISBI%nStage=13
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.6160178650170565d0
      A(3)=-0.4449487060774118d0
      A(4)=-1.0952033345276178d0
      A(5)=-1.2256030785959187d0
      A(6)=-0.2740182222332805d0
      A(7)=-0.0411952089052647d0
      A(8)=-0.1797084899153560d0
      A(9)=-1.1771530652064288d0
      A(10)=-0.4078831463120878d0
      A(11)=-0.8295636426191777d0
      A(12)=-4.7895970584252288d0
      A(13)=-0.6606671432964504d0
      B(1)=0.0271990297818803d0
      B(2)=0.1772488819905108d0
      B(3)=0.0378528418949694d0
      B(4)=0.6086431830142991d0
      B(5)=0.2154313974316100d0
      B(6)=0.2066152563885843d0
      B(7)=0.0415864076069797d0
      B(8)=0.0219891884310925d0
      B(9)=0.9893081222650993d0
      B(10)=0.0063199019859826d0
      B(11)=0.3749640721105318d0
      B(12)=1.6080235151003195d0
      B(13)=0.0961209123818189d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('TSRKC73')
      ! 7-stage third-order TSRKC73 of
      ! Toulorge, Desmet, J. Comp. Phys. 231 (2012) 2067-2091
      WRITE(*,*) 'Method TSRKC73 '
      MISBI%nStage=7
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.8083163874983830d0
      A(3)=-1.503407858773331d0
      A(4)=-1.053064525050744d0
      A(5)=-1.463149119280508d0
      A(6)=-0.6592881281087830d0
      A(7)=-1.667891931891068d0
      B(1)=0.01197052673097840d0
      B(2)=0.8886897793820711d0
      B(3)=0.4578382089261419d0
      B(4)=0.5790045253338471d0
      B(5)=0.3160214638138484d0
      B(6)=0.2483525368264122d0
      B(7)=0.06771230959408840d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE ('TSRKC84')
      ! 8-stage fourth-order TSRKC84 of
      ! Toulorge, Desmet, J. Comp. Phys. 231 (2012) 2067-2091
      WRITE(*,*) 'Method TSRKC84 '
      MISBI%nStage=8
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.7212962482279240d0
      A(3)=-0.01077336571612980d0
      A(4)=-0.5162584698930970d0
      A(5)=-1.730100286632201d0
      A(6)=-5.200129304403076d0
      A(7)=0.7837058945416420d0
      A(8)=-0.5445836094332190d0
      B(1)=0.2165936736758085d0
      B(2)=0.1773950826411583d0
      B(3)=0.01802538611623290d0
      B(4)=0.08473476372541490d0
      B(5)=0.8129106974622483d0
      B(6)=1.903416030422760d0
      B(7)=0.1314841743399048d0
      B(8)=0.2082583170674149d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE ('TSRKF84')
      ! 8-stage fourth-order TSRKF84 of
      ! Toulorge, Desmet, J. Comp. Phys. 231 (2012) 2067-2091
      WRITE(*,*) 'Method TSRKF84 '
      MISBI%nStage=8
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.5534431294501569d0 
      A(3)=0.01065987570203490d0 
      A(4)=-0.5515812888932000d0 
      A(5)=-1.885790377558741d0 
      A(6)=-5.701295742793264d0 
      A(7)=2.113903965664793d0 
      A(8)=-0.5339578826675280d0
      B(1)=0.08037936882736950d0
      B(2)=0.5388497458569843d0
      B(3)=0.01974974409031960d0
      B(4)=0.09911841297339970d0
      B(5)=0.7466920411064123d0
      B(6)=1.679584245618894d0
      B(7)=0.2433728067008188d0
      B(8)=0.1422730459001373d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('SHRK64')
      ! 6-stage fourth-order SHRK64 of
      ! Stanescu, Habashi, J. Comp. Phys. 143 (1998) 674-681
      WRITE(*,*) 'Method SHRK64'
      MISBI%nStage=6
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.4919575d0
      A(3)=-0.8946264d0
      A(4)=-1.5526678d0
      A(5)=-3.4077973d0
      A(6)=-1.0742640d0
      B(1)=0.1453095d0
      B(2)=0.4653797d0
      B(3)=0.4675397d0
      B(4)=0.7795279d0
      B(5)=0.3574327d0
      B(6)=0.15d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE ('HALERK74')
     ! 7-stage fourth-order HALERK74 of
      ! Allampalli, Hixon, Nallasamy, Sawyer, J. Comp. Phys. 228 (2009) 3837-3850
      WRITE(*,*) 'Method HALERK74'
      MISBI%nStage=7
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.647900745934d0
      A(3)=-2.704760863204d0
      A(4)=-0.460080550118d0
      A(5)=-0.500581787785d0
      A(6)=-1.906532255913d0
      A(7)=-1.450000000000d0
      B(1)=0.117322146869d0
      B(2)=0.503270262127d0
      B(3)=0.233663281658d0
      B(4)=0.283419634625d0
      B(5)=0.540367414023d0
      B(6)=0.371499414620d0
      B(7)=0.136670099385d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('HALERK64')
      ! 6-stage fourth-order HALERK64 of
      ! Allampalli, Hixon, Nallasamy, Sawyer, J. Comp. Phys. 228 (2009) 3837-3850
      MISBI%nStage=6
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A(2)=-0.691750960670d0
      A(3)=-1.727127405211d0
      A(4)=-0.694890150986d0
      A(5)=-1.039942756197d0
      A(6)=-1.531977447611d0
      B(1)=0.122000000000d0
      B(2)=0.477263056358d0
      B(3)=0.381941220320d0
      B(4)=0.447757195744d0
      B(5)=0.498614246822d0
      B(6)=0.186648570846d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE ('YRK1351')
      ! 13-stage fifth-order YRK135 of
      ! Yan, Chin. J. Chem. Phys 30 (2017) 277-286
      WRITE(*,*) 'Method YRK1351'
      MISBI%nStage=13
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      ALLOCATE(A(MISBI%nStage))
      ALLOCATE(B(MISBI%nStage))
      A( 2)=-0.33672143119427413d0
      A( 3)=-1.2018205782908164d0
      A( 4)=-2.6261919625495068d0
      A( 5)=-1.5418507843260567d0
      A( 6)=-0.2845614242371758d0
      A( 7)=-0.1700096844304301d0
      A( 8)=-1.0839412680446804d0
      A( 9)=-11.61787957751822d0
      A(10)=-4.5205208057464192d0
      A(11)=-35.86177355832474d0
      A(12)=-0.000021340899996007288d0
      A(13)=-0.066311516687861348d0
      B( 1)=0.069632640247059393d0
      B( 2)=0.088918462778092020d0
      B( 3)=1.0461490123426779d0
      B( 4)=0.42761794305080487d0
      B( 5)=0.20975844551667144d0
      B( 6)=-0.11457151862012136d0
      B( 7)=-0.01392019988507068d0
      B( 8)=4.0330655626956709d0
      B( 9)=0.35106846752457162d0
      B(10)=-0.16066651367556576d0
      B(11)=-0.0058633163225038929d0
      B(12)=0.077296133865151863d0
      B(13)=0.054301254676908338d0
      CALL LowStorageToMISBI(A,B,MISBI)
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('MRI-GARK-ERK45a')
    ! Sandu
      WRITE(*,*) 'Method: MRI-GARK-ERK45a'
      MISBI%nStage=5
      MISBI%nPhi=2
      CALL AllocateMISBI

      MISBI%a(2,1)%a(1)=1.0d0/5.0d0

      MISBI%a(3,1)%a(1)=-53.0d0/16.0d0
      MISBI%a(3,2)%a(1)=281.0d0/80.0d0

      MISBI%a(4,1)%a(1)=-36562993.0d0/71394880.0d0
      MISBI%a(4,2)%a(1)=34903117.0d0/17848720.0d0
      MISBI%a(4,3)%a(1)=-88770499.0d0/71394880.0d0

      MISBI%a(5,1)%a(1)=-7631593.0d0/71394880.0d0
      MISBI%a(5,2)%a(1)=-166232021.0d0/35697440.0d0
      MISBI%a(5,3)%a(1)=6068517.0d0/1519040.0d0
      MISBI%a(5,4)%a(1)=8644289.0d0/8924360.0d0

      MISBI%a(6,1)%a(1)=277061.0d0/303808.0d0
      MISBI%a(6,2)%a(1)=-209323.0d0/1139280.0d0
      MISBI%a(6,3)%a(1)=-1360217.0d0/1139280.0d0
      MISBI%a(6,4)%a(1)=-148789.0d0/56964.0d0
      MISBI%a(6,5)%a(1)=147889.0d0/45120.0d0


      MISBI%a(3,1)%a(2)=503.0d0/80.0d0
      MISBI%a(3,2)%a(2)=-503.0d0/80.0d0

      MISBI%a(4,1)%a(2)=-1365537.0d0/35697440.0d0
      MISBI%a(4,2)%a(2)=4963773.0d0/7139488.0d0
      MISBI%a(4,3)%a(2)=-1465833.0d0/2231090.0d0

      MISBI%a(5,1)%a(2)=66974357.0d0/35697440.0d0
      MISBI%a(5,2)%a(2)=21445367.0d0/7139488.0d0
      MISBI%a(5,3)%a(2)=-3.0d0
      MISBI%a(5,4)%a(2)=-8388609.0d0/4462180.0d0

      MISBI%a(6,1)%a(2)=-18227.0d0/7520.0d0
      MISBI%a(6,2)%a(2)=2.0d0
      MISBI%a(6,3)%a(2)=1.0d0
      MISBI%a(6,4)%a(2)=5.0d0
      MISBI%a(6,5)%a(2)=-41933.0d0/7520.0d0

      MISBI%c=0.0d0
      MISBI%c(2)=1.0d0/5.0d0
      MISBI%c(3)=2.0d0/5.0d0
      MISBI%c(4)=3.0d0/5.0d0
      MISBI%c(5)=4.0d0/5.0d0
      MISBI%c(6)=1.0d0
      MISBI%d=0.0d0
!     MISBI%d(2,1)=1.0d0
      MISBI%d(2,1)=0.0d0
      MISBI%d(3,2)=1.0d0
      MISBI%d(4,3)=1.0d0
      MISBI%d(5,4)=1.0d0
      MISBI%d(6,5)=1.0d0
      MISBI%g=0.0d0
      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0/5.0d0
      MISBI%dt(3)=1.0d0/5.0d0
      MISBI%dt(4)=1.0d0/5.0d0
      MISBI%dt(5)=1.0d0/5.0d0
      MISBI%dt(6)=1.0d0/5.0d0
    CASE('OwrenRK4')
      MISBI%nStage=5
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(2,1)%a(1)=1.0d0/2.0d0

      MISBI%a(3,2)%a(1)=1.0d0/2.0d0

      MISBI%a(4,1)%a(1)=-1.0d0/2.0d0
      MISBI%a(4,3)%a(1)=1.0d0

      MISBI%a(5,1)%a(1)=1.0d0/4.0d0
      MISBI%a(5,2)%a(1)=1.0d0/6.0d0
      MISBI%a(5,3)%a(1)=1.0d0/6.0d0
      MISBI%a(5,4)%a(1)=-1.0d0/12.0d0

      MISBI%a(6,1)%a(1)=-1.0d0/12.0d0
      MISBI%a(6,2)%a(1)=1.0d0/6.0d0
      MISBI%a(6,3)%a(1)=1.0d0/6.0d0
      MISBI%a(6,4)%a(1)=1.0d0/4.0d0

      MISBI%d(4,2)=1.0d0
      MISBI%d(6,5)=1.0d0

      MISBI%dt(1)=0.0d0
      MISBI%dt(2)=1.0d0/2.0d0
      MISBI%dt(3)=1.0d0/2.0d0
      MISBI%dt(4)=1.0d0/2.0d0
      MISBI%dt(5)=1.0d0/2.0d0
      MISBI%dt(6)=1.0d0/2.0d0

    CASE('MISBI4')
      MISBI%nStage=4
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1
      MISBI%a(2,1)%a(1)=  0.13629647842266179     
      MISBI%a(3,1)%a(1)=  0.28046239897933556     
      MISBI%a(3,2)%a(1)= -1.60351333596248577E-002
      MISBI%a(4,1)%a(1)=  0.90471335520843155     
      MISBI%a(4,2)%a(1)=  -1.0401118315403626     
      MISBI%a(4,3)%a(1)=  0.65233756348866223     
      MISBI%a(5,1)%a(1)=  6.71969845545695721E-002
      MISBI%a(5,2)%a(1)= -0.36562186260961194     
      MISBI%a(5,3)%a(1)= -0.15486147083521187     
      MISBI%a(5,4)%a(1)=  0.97036244446880304     
      MISBI%c=0.0d0
      MISBI%c(2)=  0.13629647842266179     
      MISBI%c(3)=  0.48155366095621122     
      MISBI%c(4)=  0.58375197724891603     
      MISBI%c(5)=  0.99999999992067423     
      MISBI%d=0.0d0
      MISBI%d(3,2)=  0.91409281030389122     
      MISBI%d(4,2)=   1.1427441739726025     
      MISBI%d(4,3)= -0.29521124618847322     
      MISBI%d(5,2)=  0.11296528223065300     
      MISBI%d(5,3)=  0.33736941129609332     
      MISBI%d(5,4)=  0.50374718311858491     
      MISBI%g=0.0d0
      MISBI%g(3,2)=  0.67895198329071216     
      MISBI%g(4,2)=  -1.3897416407020489     
      MISBI%g(4,3)=  0.50386457630169545     
      MISBI%g(5,2)= -0.37532860828160519     
      MISBI%g(5,3)=  0.32092502110855337     
      MISBI%g(5,4)= -0.15825968894504380     
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('MISBIC4')  
      WRITE(*,*) 'MISBIC4'
      MISBI%nStage=5
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%A(2,1)%a(1)=0.5d0

      MISBI%A(3,1)%a(1)=0.0d0
      MISBI%A(3,2)%a(1)=0.5d0

      MISBI%A(4,1)%a(1)=-0.5d0
      MISBI%A(4,2)%a(1)=0.0d0
      MISBI%A(4,3)%a(1)=1.0d0

      MISBI%A(5,1)%a(1)=0.25d0
      MISBI%A(5,2)%a(1)=1.0d0/6.0d0
      MISBI%A(5,3)%a(1)=1.0d0/6.0d0
      MISBI%A(5,4)%a(1)=-1.0d0/12.0d0

      MISBI%A(6,1)%a(1)=-1.0d0/12.0d0
      MISBI%A(6,2)%a(1)=1.0d0/6.0d0
      MISBI%A(6,3)%a(1)=1.0d0/6.0d0
      MISBI%A(6,4)%a(1)=0.25d0
      MISBI%A(6,5)%a(1)=0.0d0

      MISBI%c=0.0d0
      MISBI%d=0.0d0
      MISBI%d(4,2)=1.0d0
      MISBI%d(6,5)=1.0d0
      MISBI%g=0.0d0

      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE ('MISBIWS3')  
      MISBI%nStage=3
      MISBI%nPhi=3
      CALL AllocateMISBI
      MISBI%sPhi=0
      MISBI%A(2,1)%a(1)=0.5d0
      MISBI%sPhi(2)=1

      MISBI%A(3,1)%a(1)=-40.0d0 
      MISBI%A(3,1)%a(2)=78.0d0
      MISBI%A(3,2)%a(1)=41.0d0
      MISBI%A(3,2)%a(2)=-78.0d0
      MISBI%sPhi(3)=2

      MISBI%A(4,1)%a(1)=1.0d0/3.0d0
      MISBI%A(4,1)%a(2)=-1.0d0
      MISBI%A(4,1)%a(3)=1.0d0
      MISBI%A(4,2)%a(1)=4.0d0/3.0d0
      MISBI%A(4,2)%a(3)=-2.0d0
      MISBI%A(4,3)%a(1)=-2.0d0/3.0d0
      MISBI%A(4,3)%a(2)=1.0d0
      MISBI%A(4,3)%a(3)=1.0d0
      MISBI%sPhi(4)=3

      MISBI%c(1)=0.0d0
      MISBI%c(2)=0.5d0
      MISBI%c(3)=1.0d0
      MISBI%c(4)=1.0d0

      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
    CASE('MISTvdA') 
      MISBI%nStage=3
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1

      MISBI%A(2,1)%a(1)=  0.66666666666666696d0
      MISBI%A(3,1)%a(1)= -0.28247174703488398d0
      MISBI%A(3,2)%a(1)=  0.44444444444444398d0
      MISBI%A(4,1)%a(1)= -0.31198081960042401d0
      MISBI%A(4,2)%a(1)=  0.18082737579913699d0
      MISBI%A(4,3)%a(1)=  0.56250000000000000d0

      MISBI%c(2)=  0.66666666666666696d0     
      MISBI%c(3)=  0.66666666666666685d0     
      MISBI%c(4)=   1.0000000000000009d0     

      MISBI%d(3,2)=  0.1946360605647457d0
      MISBI%d(4,2)=  0.3971200136786614d0
      MISBI%d(4,3)=  0.2609434606211801d0

      MISBI%g(3,2)=  0.5624048933209129d0
      MISBI%g(4,2)=  0.4408467475713277d0
      MISBI%g(4,3)= -0.2459300561692391d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%A(i,j)%a(1)
        END DO
      END DO
    CASE('MISTvdB') 
      MISBI%nStage=3
      MISBI%nPhi=1
      CALL AllocateMISBI
      MISBI%sPhi=1

      MISBI%A(2,1)%a(1)=  0.66666666666666696     
      MISBI%A(3,1)%a(1)= -0.25492859100078202     
      MISBI%A(3,2)%a(1)=  0.44444444444444398     
      MISBI%A(4,1)%a(1)= -0.26452517179288798     
      MISBI%A(4,2)%a(1)=  0.11424084424766399     
      MISBI%A(4,3)%a(1)=  0.56250000000000000     
      MISBI%c(2)=  0.66666666666666696     
      MISBI%c(3)=  0.66666666666666685     
      MISBI%c(4)=   1.0000000000000009     
      MISBI%d(3,2)=  0.42668232863311001     
      MISBI%d(4,2)=  0.26570779016173801     
      MISBI%d(4,3)=  0.41489966891866698     
      MISBI%g(3,2)=  0.28904389120139701     
      MISBI%g(4,2)=  0.45113560071334202     
      MISBI%g(4,3)= -0.25006656847591002     
      MISBI%dt=0.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO

    CASE('MISBIEB4')
      WRITE(*,*) 'MISBIEB4'
      MISBI%nStage=4
      MISBI%nPhi=3
      CALL AllocateMISBI

      MISBI%A(2,1)%a(1)=0.5d0

      MISBI%A(3,1)%a(1)= 0.5d0
      MISBI%A(3,1)%a(2)=-1.0d0
      MISBI%A(3,2)%a(2)= 1.0d0

      MISBI%A(4,1)%a(1)= 1.0d0
      MISBI%A(4,1)%a(2)=-2.0d0
      MISBI%A(4,3)%a(2)= 2.0d0

      MISBI%A(5,1)%a(1)= 1.0d0
      MISBI%A(5,1)%a(2)=-3.0d0
      MISBI%A(5,1)%a(3)= 4.0d0
      MISBI%A(5,2)%a(1)= 0.0d0
      MISBI%A(5,2)%a(2)= 2.0d0
      MISBI%A(5,2)%a(3)=-4.0d0
      MISBI%A(5,3)%a(1)= 0.0d0
      MISBI%A(5,3)%a(2)= 2.0d0
      MISBI%A(5,3)%a(3)=-4.0d0
      MISBI%A(5,4)%a(1)= 0.0d0
      MISBI%A(5,4)%a(2)=-1.0d0
      MISBI%A(5,4)%a(3)= 4.0d0

      MISBI%c(1)=0.0d0
      MISBI%c(2)=0.5d0
      MISBI%c(3)=0.5d0
      MISBI%c(4)=1.0d0
      MISBI%c(5)=1.0d0
      DO i=2,MISBI%nStage+1
        DO j=1,i-1
          MISBI%dt(i)=MISBI%dt(i)+MISBI%a(i,j)%a(1)
        END DO
      END DO
      DO i=2,MISBI%nStage+1
        Fac=1.0d0
        DO k=1,MISBI%nPhi
          DO j=1,i-1
            MISBI%a(i,j)%a(k)=MISBI%a(i,j)%a(k)*Fac
          END DO
!         Fac=Fac*MISBI%dt(i)/FLOAT(k)
          Fac=Fac/FLOAT(k)
        END DO
      END DO
  END SELECT    
CONTAINS  
SUBROUTINE LowStorageToMISBI(A,B,MISBI)
  REAL(8) :: A(:)
  REAL(8) :: B(:)
  TYPE(MISBIMethod_T) :: MISBI

  INTEGER :: i,j

  DO i=1,MISBI%nStage
    MISBI%d(i+1,i)=1.0d0
  END DO  
  DO i=1,MISBI%nStage
    MISBI%A(i+1,i)%a(1)=B(i)
    DO j=i-1,1,-1
      MISBI%A(i+1,j)%a(1)=A(j+1)*MISBI%A(i+1,j+1)%a(1)+B(j)
    END DO
  END DO
  DO i=MISBI%nStage,1,-1
    DO j=1,MISBI%nStage
      MISBI%A(i+1,j)%a=MISBI%A(i+1,j)%a-MISBI%A(i,j)%a
    END DO  
  END DO  
END SUBROUTINE LowStorageToMISBI

SUBROUTINE AllocateMISBI
  ALLOCATE(MISBI%c(MISBI%nStage+1))
  MISBI%c=0.0d0
  ALLOCATE(MISBI%d(MISBI%nStage+1,MISBI%nStage))
  MISBI%d=0.0d0
  ALLOCATE(MISBI%g(MISBI%nStage+1,MISBI%nStage))
  MISBI%g=0.0d0
  ALLOCATE(MISBI%dt(MISBI%nStage+1))
  MISBI%dt=0.0d0
  ALLOCATE(MISBI%sPhi(MISBI%nStage+1))
  ALLOCATE(MISBI%a(MISBI%nStage+1,MISBI%nStage))
  DO i=1,MISBI%nStage+1
    DO j=1,MISBI%nStage
      ALLOCATE(MISBI%A(i,j)%a(MISBI%nPhi))
      MISBI%A(i,j)%a=0.0d0
    END DO 
  END DO 
END SUBROUTINE AllocateMISBI

END SUBROUTINE MISBIMethod

END MODULE MISBI_Mod
