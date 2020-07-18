MODULE BSeries_Mod
  USE BiColourTree_Mod
  USE PowerSeries_Mod
  USE Polynom_Mod
  IMPLICIT NONE

  TYPE Series_T
    REAL(8), ALLOCATABLE :: a(:)
  END TYPE Series_T

  TYPE PSeries_T
    TYPE(PKoeff_T), ALLOCATABLE :: a(:)
  END TYPE PSeries_T

  TYPE BSeries_T
    TYPE(Series_T), POINTER :: Phi(:)=>NULL()
  END TYPE BSeries_T

  TYPE BPSeries_T
    TYPE(PSeries_T), POINTER :: Phi(:)=>NULL()
  END TYPE BPSeries_T

  TYPE(BSeries_T) :: E
  TYPE(BSeries_T) :: D
  TYPE(BSeries_T) :: Db
  TYPE(BSeries_T) :: Dw
CONTAINS

SUBROUTINE MultPowerSeries(PB,Pow,B,Order)
  TYPE(BSeries_T) :: PB
  TYPE(PowerSeries_T) :: Pow
  TYPE(BSeries_T) :: B
  INTEGER :: Order

  INTEGER :: i,iLoc,j,p,pLoc
  TYPE(TreeP_T), POINTER :: T
  REAL(8) :: a(SIZE(PB%Phi(Order)%a))

  a=0.0d0
  DO p=0,Order
    T=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      DO j=0,T%TP%NumberPowers
        pLoc=T%TP%Powers(j)%TP%Order
        IF (pLoc==Order) THEN
          iLoc=T%TP%Powers(j)%TP%NumberOrder
          a(iLoc)=a(iLoc)+B%Phi(p)%a(i)*Pow%Koeff(j)
        END IF
      END DO
      T=>T%Next
    END DO
  END DO
  PB%Phi(Order)%a=a

END SUBROUTINE MultPowerSeries
  
SUBROUTINE MultPowerSeries1(PB,Pow,B,Order)
  TYPE(BSeries_T) :: PB
  TYPE(PowerSeries_T) :: Pow
  TYPE(BSeries_T) :: B
  INTEGER :: Order

  INTEGER :: i,iLoc,j,p,pLoc
  TYPE(TreeP_T), POINTER :: T

  DO p=1,Order
    PB%Phi(p)%a=0.0d0
  END DO  
  DO p=0,Order
    T=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      DO j=0,T%TP%NumberPowers
        pLoc=T%TP%Powers(j)%TP%Order
        IF (pLoc<=Order) THEN
          iLoc=T%TP%Powers(j)%TP%NumberOrder
          PB%Phi(pLoc)%a(iLoc)=PB%Phi(pLoc)%a(iLoc)+B%Phi(p)%a(i)*Pow%Koeff(j)
        END IF
      END DO
      T=>T%Next
    END DO
  END DO

END SUBROUTINE MultPowerSeries1

SUBROUTINE AllocateBSeries(B,pMax)
  INTEGER :: pMax
  TYPE(BSeries_T) :: B

  INTEGER :: p
  IF (.NOT.ASSOCIATED(B%Phi)) THEN
    ALLOCATE(B%Phi(0:pMax))
    DO p=0,pMax
      ALLOCATE(B%Phi(p)%a(ListTree(p)%LenListTree))
      B%Phi(p)%a=0.0d0
    END DO  
  END IF
END SUBROUTINE AllocateBSeries

SUBROUTINE AllocateBPSeries(B,pMax)
  INTEGER :: pMax
  TYPE(BPSeries_T) :: B

  INTEGER :: i,p
  IF (.NOT.ASSOCIATED(B%Phi)) THEN
    ALLOCATE(B%Phi(0:pMax))
    DO p=0,pMax
      ALLOCATE(B%Phi(p)%a(ListTree(p)%LenListTree))
      DO i=1,ListTree(p)%LenListTree
        ALLOCATE(B%Phi(p)%a(i)%Koeff(0:0))
        B%Phi(p)%a(i)%Grad=0
        B%Phi(p)%a(i)%Koeff=0.0d0
      END DO  
    END DO  
  END IF
END SUBROUTINE AllocateBPSeries

SUBROUTINE InitBSeries(pMax)
  INTEGER :: pMax

  INTEGER :: i,p
  REAL(8) :: Gam
  TYPE(TreeP_T), POINTER :: T1

  ALLOCATE(E%Phi(0:pMax))
  ALLOCATE(D%Phi(0:pMax))
  ALLOCATE(Db%Phi(0:pMax))
  ALLOCATE(Dw%Phi(0:pMax))
  DO p=0,pMax
    ALLOCATE(E%Phi(p)%a(ListTree(p)%LenListTree))
    ALLOCATE(D%Phi(p)%a(ListTree(p)%LenListTree))
    ALLOCATE(Db%Phi(p)%a(ListTree(p)%LenListTree))
    ALLOCATE(Dw%Phi(p)%a(ListTree(p)%LenListTree))
  END DO  
  E%Phi(0)%a(1)=1.0d0
  D%Phi(0)%a(1)=0.0d0
  Db%Phi(0)%a(1)=0.0d0
  Dw%Phi(0)%a(1)=0.0d0
  DO p=1,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      Gam=Gamma(T1%TP)
      E%Phi(p)%a(i)=1.0d0/Gam
      D%Phi(p)%a(i)=0.0d0
      Db%Phi(p)%a(i)=0.0d0
      Dw%Phi(p)%a(i)=0.0d0
      T1=>T1%Next
    END DO  
  END DO  
  D%Phi(1)%a(1)=1.0d0
  Db%Phi(1)%a(1)=1.0d0
  Dw%Phi(1)%a(2)=1.0d0

END SUBROUTINE InitBSeries

SUBROUTINE ComposeBSeries(AB,A,B,pMax)
  TYPE(BSeries_T) :: AB,A,B
  INTEGER :: pMax

  INTEGER :: i,k,l,p
  INTEGER :: iLocA,pLocA
  INTEGER :: iLocB,pLocB
  REAL(8) :: Temp
  TYPE(TreeP_T), POINTER :: T1
  TYPE(TreeP_T), POINTER :: T2
  TYPE(STreeP_T), POINTER :: STree

  CALL AllocateBSeries(AB,pMax)
  DO p=0,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      AB%Phi(p)%a(i)=0.0d0
      STree=>T1%TP%STree
      k=1
      DO 
        IF (ASSOCIATED(STree)) THEN
          pLocB=STree%TP%Order
          iLocB=STree%TP%NumberOrder
          Temp=1.0d0
          DO l=1,SIZE(STree%ListCuts)
            pLocA=NumberToTree(STree%ListCuts(l))%TP%Order
            iLocA=NumberToTree(STree%ListCuts(l))%TP%NumberOrder
            Temp=Temp*A%Phi(pLocA)%a(iLocA)**STree%NumberListCuts(l) 
          END DO  
          AB%Phi(p)%a(i)=AB%Phi(p)%a(i)+B%Phi(pLocB)%a(iLocB)*Temp
          k=k+1
          STree=>Stree%Next
        ELSE  
          EXIT
        END IF
      END DO
      T1=>T1%Next
    END DO  
  END DO  
END SUBROUTINE ComposeBSeries

END MODULE BSeries_Mod
