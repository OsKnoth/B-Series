MODULE Polynom_Mod

  IMPLICIT NONE

  TYPE PKoeff_T
    INTEGER :: Grad
    REAL(8), ALLOCATABLE :: Koeff(:)
  END TYPE PKoeff_T

CONTAINS  

SUBROUTINE CopyPolynom(P1,P2)
  TYPE(PKoeff_T), INTENT(OUT) :: P1
  TYPE(PKoeff_T), INTENT(IN) :: P2

  IF (ALLOCATED(P1%Koeff)) THEN
    DEALLOCATE(P1%Koeff)
  END IF
  P1%Grad=P2%Grad
  ALLOCATE(P1%Koeff(0:P1%Grad))
  P1%Koeff=P2%Koeff
END SUBROUTINE CopyPolynom

SUBROUTINE MultPolynom(P,P1,P2)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  TYPE(PKoeff_T), INTENT(IN) :: P1,P2

  INTEGER :: i,s,GradP,GradP1,GradP2
  REAL(8) :: KoeffP1(0:P1%Grad+P2%Grad)
  REAL(8) :: KoeffP2(0:P1%Grad+P2%Grad)

  GradP1=P1%Grad
  GradP2=P2%Grad
  KoeffP1=0.0d0
  KoeffP2=0.0d0
  KoeffP1(0:P1%Grad)=P1%Koeff(0:GradP1)
  KoeffP2(0:P2%Grad)=P2%Koeff(0:GradP2)
  GradP=GradP1+GradP2

  IF (ALLOCATED(P%Koeff)) THEN
    IF (GradP>SIZE(P%Koeff(1:))) THEN
      DEALLOCATE(P%Koeff)
      ALLOCATE(P%Koeff(0:GradP))
    END IF
  ELSE
    ALLOCATE(P%Koeff(0:GradP))
  END IF
  P%Koeff=0.0d0
  DO s=0,GradP
    DO i=0,s
      P%Koeff(s)=P%Koeff(s)+KoeffP1(i)*KoeffP2(s-i)
    END DO
  END DO  
  P%Grad=GradP
END SUBROUTINE MultPolynom

SUBROUTINE AddPolynom(P,P1,P2)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  TYPE(PKoeff_T), INTENT(IN) :: P1,P2

  INTEGER :: s,GradP,GradP1,GradP2
  REAL(8) :: KoeffP1(0:P1%Grad+P2%Grad)
  REAL(8) :: KoeffP2(0:P1%Grad+P2%Grad)

  GradP1=P1%Grad
  GradP2=P2%Grad
  KoeffP1=0.0d0
  KoeffP2=0.0d0
  KoeffP1(0:GradP1)=P1%Koeff(0:GradP1)
  KoeffP2(0:GradP2)=P2%Koeff(0:GradP2)
  GradP=MAX(GradP1,GradP2)

  IF (ALLOCATED(P%Koeff)) THEN
    IF (GradP>SIZE(P%Koeff(1:))) THEN
      DEALLOCATE(P%Koeff)
      ALLOCATE(P%Koeff(0:GradP))
    END IF
  ELSE
    ALLOCATE(P%Koeff(0:GradP))
  END IF
  P%Koeff=0.0d0
  DO s=0,GradP
    P%Koeff(s)=KoeffP1(s)+KoeffP2(s)
  END DO  
  P%Grad=GradP
END SUBROUTINE AddPolynom

SUBROUTINE AddPolynomScalar(P,P1,Val,P2)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  TYPE(PKoeff_T), INTENT(IN) :: P1
  REAL(8) :: Val
  TYPE(PKoeff_T), INTENT(IN) :: P2

  INTEGER :: s,GradP,GradP1,GradP2
  REAL(8) :: KoeffP1(0:P1%Grad+P2%Grad)
  REAL(8) :: KoeffP2(0:P1%Grad+P2%Grad)

  GradP1=P1%Grad
  GradP2=P2%Grad
  KoeffP1=0.0d0
  KoeffP2=0.0d0
  KoeffP1(0:GradP1)=P1%Koeff(0:GradP1)
  KoeffP2(0:GradP2)=P2%Koeff(0:GradP2)
  GradP=MAX(GradP1,GradP2)

  IF (ALLOCATED(P%Koeff)) THEN
    IF (GradP>SIZE(P%Koeff(1:))) THEN
      DEALLOCATE(P%Koeff)
      ALLOCATE(P%Koeff(0:GradP))
    END IF
  ELSE
    ALLOCATE(P%Koeff(0:GradP))
  END IF
  P%Koeff=0.0d0
  DO s=0,GradP
    P%Koeff(s)=KoeffP1(s)+Val*KoeffP2(s)
  END DO  
  P%Grad=GradP
END SUBROUTINE AddPolynomScalar

SUBROUTINE AddScalarPolynom(P,Val,Pos)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  REAL(8) :: Val
  INTEGER :: Pos

  INTEGER :: GradP
  REAL(8) :: PKoeff(0:MAX(Pos,P%Grad))

  GradP=MAX(Pos,P%Grad)
  IF (ALLOCATED(P%Koeff)) THEN
    IF (GradP>SIZE(P%Koeff(1:))) THEN
      PKoeff=0.0d0
      PKoeff(0:P%Grad)=P%Koeff(0:P%Grad)
      DEALLOCATE(P%Koeff)
      ALLOCATE(P%Koeff(0:GradP))
      P%Grad=GradP
      P%Koeff(0:P%Grad)=PKoeff(0:P%Grad)
    END IF
  ELSE
    ALLOCATE(P%Koeff(0:GradP))
    P%Koeff=0.0d0
  END IF

  P%Koeff(Pos)=P%Koeff(Pos)+Val
  P%Grad=GradP

END SUBROUTINE AddScalarPolynom   

SUBROUTINE MultScalarPolynom(P,Val)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  REAL(8) :: Val

  IF (ALLOCATED(P%Koeff)) THEN
    P%Koeff(0:P%Grad)=P%Koeff(0:P%Grad)*Val
  ELSE
    ALLOCATE(P%Koeff(0:0))
    P%Grad=0
    P%Koeff=0.0d0
  END IF

END SUBROUTINE MultScalarPolynom   

SUBROUTINE IntegratePolynom(P,P1)
  TYPE(PKoeff_T), INTENT(INOUT) :: P
  TYPE(PKoeff_T), INTENT(IN) :: P1

  INTEGER :: s,GradP
  REAL(8) :: KoeffP1(0:P1%Grad)


  KoeffP1=P1%Koeff
  GradP=P1%Grad+1

  IF (ALLOCATED(P%Koeff)) THEN
    IF (GradP>SIZE(P%Koeff(1:))) THEN
      DEALLOCATE(P%Koeff)
      ALLOCATE(P%Koeff(0:GradP))
    END IF
  ELSE
    ALLOCATE(P%Koeff(0:GradP))
  END IF
  P%Koeff=0.0d0
  DO s=1,GradP
    P%Koeff(s)=KoeffP1(s-1)/FLOAT(s)
  END DO  
  P%Grad=GradP

END SUBROUTINE IntegratePolynom

SUBROUTINE UnitPolynom(P)
  TYPE(PKoeff_T), INTENT(INOUT) :: P

  IF (ALLOCATED(P%Koeff)) THEN
    DEALLOCATE(P%Koeff)
  END IF
  P%Grad=0
  ALLOCATE(P%Koeff(0:0))
  P%Koeff(0)=1.0d0
END SUBROUTINE UnitPolynom

SUBROUTINE PowerOnePolynom1(P)
  TYPE(PKoeff_T), INTENT(INOUT) :: P

  IF (ALLOCATED(P%Koeff)) THEN
    DEALLOCATE(P%Koeff)
  END IF
  P%Grad=1
  ALLOCATE(P%Koeff(0:P%Grad))
  P%Koeff(0)=0.0d0
  P%Koeff(1)=1.0d0
END SUBROUTINE PowerOnePolynom1

SUBROUTINE PowerOnePolynom(P)
  TYPE(PKoeff_T), INTENT(INOUT) :: P(:)

  INTEGER :: i

  DO i=1,UBOUND(P,1)
    IF (ALLOCATED(P(i)%Koeff)) THEN
      DEALLOCATE(P(i)%Koeff)
    END IF
    P(i)%Grad=1
    ALLOCATE(P(i)%Koeff(0:P(i)%Grad))
    P(i)%Koeff(0)=0.0d0
    P(i)%Koeff(1)=1.0d0/FLOAT(i)
  END DO
END SUBROUTINE PowerOnePolynom

FUNCTION ValuePolynom(P,dt)
  REAL(8) :: ValuePolynom
  TYPE(PKoeff_T) :: P
  REAL(8) :: dt

  INTEGER :: s
  ValuePolynom=P%Koeff(P%Grad)
  DO s=P%Grad-1,0,-1 
    ValuePolynom=ValuePolynom*dt+P%Koeff(s)
  END DO  
END FUNCTION ValuePolynom

END MODULE Polynom_Mod
