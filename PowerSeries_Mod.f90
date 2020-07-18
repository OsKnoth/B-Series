MODULE PowerSeries_Mod

  IMPLICIT NONE

  TYPE PowerSeries_T
    REAL(8), POINTER :: Koeff(:)=>NULL()
  END TYPE PowerSeries_T

  TYPE(PowerSeries_T) :: Inv
  TYPE(PowerSeries_T) :: Expo
  TYPE(PowerSeries_T), ALLOCATABLE :: Phi(:)
  INTEGER :: PhiMax=4

CONTAINS

SUBROUTINE InvPow(Pow,Gamma,Order)
  TYPE(PowerSeries_T) :: Pow
  REAL(8) :: Gamma
  INTEGER :: Order

  INTEGER :: i

  IF (.NOT.ASSOCIATED(Pow%Koeff)) THEN
    ALLOCATE(Pow%Koeff(0:Order))
  END IF  
  Pow%Koeff(0)=1.0d0
  DO i=1,Order
    Pow%Koeff(i)=Pow%Koeff(i-1)*Gamma
  END DO  
  
END SUBROUTINE InvPow

SUBROUTINE PhiPow(Pow,Gamma,Order,Number)
  TYPE(PowerSeries_T) :: Pow
  REAL(8) :: Gamma
  INTEGER :: Order
  INTEGER :: Number

  INTEGER :: i

  IF (.NOT.ASSOCIATED(Pow%Koeff)) THEN
    ALLOCATE(Pow%Koeff(0:Order))
  END IF  
  Pow%Koeff(0)=1.0d0
  DO i=1,Number
    Pow%Koeff(0)=Pow%Koeff(0)/FLOAT(i)
  END DO  
  DO i=1,Order
    Pow%Koeff(i)=Pow%Koeff(i-1)*Gamma/FLOAT(i+Number)
  END DO  
  
END SUBROUTINE PhiPow

SUBROUTINE JPow(Pow,Order)
  TYPE(PowerSeries_T) :: Pow
  INTEGER :: Order

  IF (.NOT.ASSOCIATED(Pow%Koeff)) THEN
    ALLOCATE(Pow%Koeff(0:Order))
  END IF  
  Pow%Koeff=0.0d0
  Pow%Koeff(1)=1.0d0
END SUBROUTINE JPow

SUBROUTINE InitPowerSeries(Order)
  INTEGER :: Order
  INTEGER :: i

  ALLOCATE(Inv%Koeff(0:Order))
  DO i=0,Order
    Inv%Koeff(i)=1.0d0
  END DO

  ALLOCATE(Expo%Koeff(0:Order))
  DO i=0,Order
    Expo%Koeff(i)=1.0d0
  END DO

END SUBROUTINE InitPowerSeries

END MODULE PowerSeries_Mod

