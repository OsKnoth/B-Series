MODULE CompareAlgorithm_Mod
  USE MISBI_Mod
  USE BiColourTree_Mod
  USE BSeries_Mod
  
  IMPLICIT NONE
CONTAINS 
SUBROUTINE OrderWenschKnoth(MISBI,Phi2MISBI,Exact,pMax)

  IMPLICIT NONE
  TYPE(MISBIMethod_T) :: MISBI
  TYPE(BSeries_T) :: Phi2MISBI(:)
  TYPE(BSeries_T) :: Exact
  INTEGER :: pMax

  INTEGER :: Stages
  REAL(8) :: alpha(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: g(MISBI%nStage+1,MISBI%nStage+1)
  TYPE(Koeff_T) :: beta(MISBI%nStage+1,MISBI%nStage)

  REAL(8) :: R(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: A(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: DD(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: gTemp(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: betaTemp(MISBI%nStage+1,MISBI%nStage+1)
  REAL(8) :: b(1,MISBI%nStage+1)
  REAL(8) :: bTilde(1,MISBI%nStage+1)
  REAL(8) :: ww(1,MISBI%nStage+1)
  REAL(8) :: w1(1,MISBI%nStage+1)
  REAL(8) :: c(MISBI%nStage+1,1)
  REAL(8) :: c2(MISBI%nStage+1,1)
  REAL(8) :: wc(MISBI%nStage+1,1)
  REAL(8) :: cTilde(MISBI%nStage+1,1)
  REAL(8) :: dt(MISBI%nStage+1)
  REAL(8) :: Ord(30)
  CHARACTER(LenTreeString) :: String
  INTEGER :: Pos(2)
  INTEGER :: pT,iT

  INTEGER :: info
  INTEGER :: i,j,k

  Stages=MISBI%nStage
  beta=MISBI%a
  betaTemp=0.0d0
  alpha=0.0d0
  g=0.0d0
  gTemp=0.0d0
  DO i=1,Stages+1
    dt(i)=MISBI%dt(i)
    DO j=1,Stages
      alpha(i,j)=MISBI%d(i,j)
      g(i,j)=MISBI%g(i,j)
      gTemp(i,j)=g(i,j)
      betaTemp(i,j)=beta(i,j)%a(1)
    END DO
  END DO  

  R=0.0d0
  DD=0.0d0
  DO i=1,Stages+1
    R(i,i)=1.0d0
    DD(i,i)=dt(i)
  END DO  
  A=0.0d0
  DO i=2,Stages+1
    A(i,i)=DD(i,i)
    DO j=2,i-1
      R(i,j)=R(i,j)-alpha(i,j)-gTemp(i,j)
    END DO  
  END DO  
  DO i=1,Stages+1
    DO k=1,Stages+1
      A(i,k)=A(i,k)/R(i,i)
    END DO
    DO j=i+1,Stages+1
      DO k=1,Stages+1
        A(j,k)=A(j,k)-R(j,i)*A(i,k)
      END DO  
    END DO  
  END DO  
  bTilde(1:1,:)=A(Stages+1:Stages+1,:)

  DO i=1,Stages+1
    DO j=1,Stages+1
      A(i,j)=betaTemp(i,j)
    END DO
  END DO
  DO i=1,Stages+1
    DO k=1,Stages+1
      A(i,k)=A(i,k)/R(i,i)
    END DO
    DO j=i+1,Stages+1
      DO k=1,Stages+1
        A(j,k)=A(j,k)-R(j,i)*A(i,k)
      END DO
    END DO  
  END DO  
  b(1:1,:)=A(Stages+1:Stages+1,:)

  DO i=1,Stages+1
    c(i,1)=0.0d0
    DO j=1,Stages+1
      c(i,1)=c(i,1)+A(i,j)
    END DO  
  END DO  
  DO i=1,Stages+1
    c2(i,1)=c(i,1)*c(i,1)
  END DO  
! cTilde=Matmul(alpha,c)
  cTilde=0.0d0
  DO i=2,Stages+1
    DO j=2,i-1
      cTilde(i,1)=cTilde(i,1)+alpha(i,j)*c(j,1)
    END DO
  END DO  
!  First Order
!     'Order 1 C'
!      Order=SUM(b)
      Ord(1)=SUM(b)
      String='.'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 1',1,Ord(1),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!
!   Second Ord  
!    'Order 2 C'
!     Ord=SUM(Matmul(b,c)-0.5d0) 
      Ord(2)=0.0d0
      DO i=1,Stages+1
        Ord(2)=Ord(2)+b(1,i)*c(i,1)
      END DO  
      String='[.]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 2',1,Ord(2),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 2 E'
!     Ord=SUM(Matmul(bTilde,(c+cTilde))-1.0d00) 
      Ord(3)=0.0d0
      DO i=1,Stages+1
        Ord(3)=Ord(3)+bTilde(1,i)*(c(i,1)+cTilde(i,1))
      END DO  
      String='[o]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 2',2,Ord(3)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

      WRITE(*,*)

!      
!  Third Order  
!
!    'Order 3 C',Ord(4)
!     Ord=SUM(Matmul(b,c2)-1.0d0/3.0d0) 
      Ord(4)=SUM(Matmul(b,c2)) 
      String='[.,.]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',1,Ord(4),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 3 C',Ord(5)
      Ord(5)=SUM(Matmul(Matmul(b,A),c))
      String='[[.]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',2,Ord(5),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 3 E',Ord(6)
      ww=Matmul(bTilde,A)
      Ord(6)=SUM(Matmul(ww,c))
      ww=Matmul(bTilde,alpha)
      ww=Matmul(ww,A)
      Ord(6)=Ord(6)+SUM(Matmul(ww,c))
      String='[(.)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',3,0.5d0*Ord(6),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!     Ord=Ord-1.0d0/3.0d0
!    'Order 3 E',Ord(7)
      ww=Matmul(bTilde,alpha+0.5d0*gTemp)
!     ww=Matmul(ww,R)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(7)=3.0d0*SUM(Matmul(ww,c+cTilde))
      ww=Matmul(bTilde,DD)
      Ord(7)=Ord(7)+SUM(Matmul(ww,c+2.0d0*cTilde))
      String='[(o)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',4,Ord(7)/6.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 3 E',Ord(8)
!     ww=Matmul(b,R)
      ww=b
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(8)=SUM(Matmul(ww,c+cTilde))
      String='[[o]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',5,0.5d0*Ord(8),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 3 E',Ord(9)
      DO i=1,Stages+1
        c2(i,1)=c2(i,1)+cTilde(i,1)**2+c(i,1)*cTilde(i,1)
      END DO  
      Ord(9)=SUM(Matmul(bTilde,c2))
      String='[o,o]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 3',6,Ord(9)/3.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

      WRITE(*,*)

!      
!  Fourth Order  
!
!    'Order 4 C',Ord(11)     1
      Ord(11)=0.0d0
      DO i=1,Stages+1
        Ord(11)=Ord(11)+b(1,i)*c(i,1)**3
      END DO  
      String='[.,.,.]'  ! f'''FFF
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',1,Ord(11),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

!    'Order 4 E',Ord(12)    2
      Ord(12)=0.0d0
      DO i=1,Stages+1
        Ord(12)=Ord(12)+bTilde(1,i)*(c(i,1)**2+cTilde(i,1)**2)*(c(i,1)+cTilde(i,1))
      END DO  
      String='[o,o,o]' ! g'''FFF
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',2,Ord(12)/4.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

!    'Order 4 C',Ord(13)     3
      Ord(13)=0.0d0
      wc=Matmul(A,c)
      DO i=1,Stages+1
        Ord(13)=Ord(13)+b(1,i)*wc(i,1)*c(i,1)
      END DO  
      String='[.,[.]]' ! f''f'FF
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',3,Ord(13),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

!    'Order 4 E',Ord(14)    4
      Ord(14)=0.0d0
      DO i=1,Stages+1
        ww(1,i)=b(1,i)*c(i,1)
      END DO  
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(14)=Ord(14)+ww(1,i)*(c(i,1)+cTilde(i,1))
      END DO  
      String='[.,[o]]'  
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',4,Ord(14)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

!    'Order 4 E',Ord(15)
      Ord(15)=0.0d0
      DO i=1,Stages+1
        ww(1,i)=bTilde(1,i)*(2.0d0*c(i,1)+cTilde(i,1))
      END DO  
      ww=Matmul(ww,A)
      DO i=1,Stages+1
        Ord(15)=Ord(15)+ww(1,i)*c(i,1)
      END DO  
      DO i=1,Stages+1
        ww(1,i)=bTilde(1,i)*(c(i,1)+2.0d0*cTilde(i,1))
      END DO  
      ww=Matmul(ww,alpha)
      ww=Matmul(ww,A)
      DO i=1,Stages+1
        Ord(15)=Ord(15)+ww(1,i)*c(i,1)
      END DO  
      String='[o,(.)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',5,Ord(15)/6.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(16)
      Ord(16)=0.0d0
      DO i=1,Stages+1
        ww(1,i)=2.0d0*bTilde(1,i)*(c(i,1)+cTilde(i,1))
      END DO  
      ww=Matmul(ww,alpha)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(16)=Ord(16)+ww(1,i)*(c(i,1)+cTilde(i,1))
      END DO  

      DO i=1,Stages+1
        ww(1,i)=4.0d0/3.0d0*bTilde(1,i)*(c(i,1)+0.5d0*cTilde(i,1))
      END DO  
      ww=Matmul(ww,gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(16)=Ord(16)+ww(1,i)*(c(i,1)+cTilde(i,1))
      END DO  

      DO i=1,Stages+1
        ww(1,i)=bTilde(1,i)*(c(i,1)+1.0d0/3.0d0*cTilde(i,1))
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(16)=Ord(16)+ww(1,i)*c(i,1)
      END DO  

      DO i=1,Stages+1
        ww(1,i)=bTilde(1,i)*(5.0d0/3.0d0*c(i,1)+cTilde(i,1))
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(16)=Ord(16)+ww(1,i)*cTilde(i,1)
      END DO  

      String='[o,(o)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',6,Ord(16)/8.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 C',Ord(17)
      Ord(17)=0.0d0
      ww=MATMUL(b,A)
      DO i=1,Stages+1
        Ord(17)=Ord(17)+ww(1,i)*c(i,1)**2
      END DO  

      String='[[.,.]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',7,Ord(17),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(18)
      Ord(18)=0.0d0
      ww=b
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      DO i=1,Stages+1
        Ord(18)=Ord(18)+ww(1,i)*(c(i,1)**2+c(i,1)*cTilde(i,1)+cTilde(i,1)**2)
      END DO  

      String='[[o,o]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',8,Ord(18)/3.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(19)
      Ord(19)=0.0d0
      ww=MATMUL(bTilde,alpha)+bTilde
      ww=MATMUL(ww,A)
      Ord(19)=SUM(MATMUL(ww,c*c))

      String='[(.,.)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',9,Ord(19)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(20)
      ww=Matmul(bTilde,alpha+0.5d0*gTemp)
!     ww=Matmul(ww,R)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(20)=4.0d0*SUM(Matmul(ww,c*c+c*cTilde+cTilde*cTilde))
      ww=Matmul(bTilde,DD)
      Ord(20)=Ord(20)+SUM(Matmul(ww,c*c+2.0d0*c*cTilde+3.0d0*cTilde*cTilde))
      String='[(o,o)]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',10,Ord(20)/12.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String) 

!    'Order 4 C',Ord(21)
      Ord(21)=0.0d0
      ww=MATMUL(b,A)
      ww=MATMUL(ww,A)
      Ord(21)=SUM(MATMUL(ww,c))

      String='[[[.]]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',11,Ord(21),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(22)
      Ord(22)=0.0d0
      ww=MATMUL(b,A)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(22)=SUM(MATMUL(ww,c+cTilde))

      String='[[[o]]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',12,Ord(22)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(23)
      Ord(23)=0.0d0
      ww=b
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=Matmul(ww,alpha)+ww
      ww=Matmul(ww,A)
      Ord(23)=SUM(MATMUL(ww,c))

      String='[[(.)]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',13,Ord(23)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(24)
      Ord(24)=0.0d0
      ww=3.0d0*b
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=Matmul(ww,alpha+0.5d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(24)=SUM(MATMUL(ww,c+cTilde))

      ww=b
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=Matmul(ww,DD)
      Ord(24)=Ord(24)+SUM(MATMUL(ww,c+2.0d0*cTilde))

      String='[[(o)]]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',14,Ord(24)/6.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(25)
      Ord(25)=0.0d0
      ww=MATMUL(bTilde,alpha)+bTilde
      ww=MATMUL(ww,A)
      ww=MATMUL(ww,A)
      Ord(25)=SUM(MATMUL(ww,c))

      String='[([.])]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',15,Ord(25)/2.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(26)
      Ord(26)=0.0d0
      ww=MATMUL(bTilde,alpha)+bTilde
      ww=MATMUL(ww,A)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(26)=SUM(MATMUL(ww,c+cTilde))

      String='[([o])]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',16,Ord(26)/4.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(27)
      Ord(27)=0.0d0
      ww=MATMUL(3.0d0*bTilde,alpha+0.5d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=MATMUL(ww,alpha)+ww
      ww=Matmul(ww,A)
      Ord(27)=SUM(MATMUL(ww,c))
      ww=MATMUL(bTilde,DD)
      ww=MATMUL(2.0d0*ww,alpha)+ww
      ww=Matmul(ww,A)
      Ord(27)=Ord(27)+SUM(MATMUL(ww,c))

      String='[((.))]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',17,Ord(27)/6.0d0,Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)

!    'Order 4 E',Ord(28)
      Ord(28)=0.0d0
      ww=MATMUL(0.5d0*bTilde,alpha+0.5d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=MATMUL(ww,alpha+0.5d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(28)=Ord(28)+SUM(MATMUL(ww,c+cTilde))

      ww=MATMUL(1.0d0/6.0d0*bTilde,alpha+0.5d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      ww=Matmul(ww,DD)
      Ord(28)=Ord(28)+SUM(MATMUL(ww,c+2.0d0*cTilde))

      ww=MATMUL(1.0d0/4.0d0*bTilde,DD)
      ww=MATMUL(ww,alpha+1.0d0/3.0d0*gTemp)
      DO i=Stages+1,1,-1
        ww(1,i)=ww(1,i)/R(i,i)
        DO j=1,i-1
          ww(1,j)=ww(1,j)-R(i,j)*ww(1,i)
        END DO  
      END DO  
      ww=Matmul(ww,DD)
      Ord(28)=Ord(28)+SUM(MATMUL(ww,c+cTilde))

      ww=MATMUL(1.0d0/24.0d0*bTilde,DD)
      ww=Matmul(ww,DD)
      Ord(28)=Ord(28)+SUM(MATMUL(ww,c+3.0d0*cTilde))

      String='[((o))]'
      Pos=FindTreeStringPos(String)
      pT=Pos(1)
      iT=Pos(2)
      WRITE(*,*) 'Order 4',18,Ord(28),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),Exact%Phi(pT)%a(iT),TRIM(String)  

END SUBROUTINE OrderWenschKnoth

SUBROUTINE OrderConditionsOwren(MISBI,Phi2MISBI,pMax)
  TYPE(MISBIMethod_T) :: MISBI
  TYPE(BSeries_T) :: Phi2MISBI(:)
  INTEGER :: pMax

  INTEGER :: i,j,k,l
  INTEGER :: nS,nPhi
  REAL(8) :: alpha(MISBI%nStage,MISBI%nStage,pMax)
  REAL(8) :: beta(MISBI%nStage,pMax)
  REAL(8) :: c(MISBI%nStage)
  REAL(8) :: kFac,lFac
  INTEGER :: sPhi(MISBI%nStage+1)
  REAL(8) :: Order(30)
  CHARACTER(LenTreeString) :: String(30)
  INTEGER :: Pos(2)
  INTEGER :: pT,iT

  REAL(8) :: TT

  nS=MISBI%nStage
  nPhi=MISBI%nStage
  sPhi=MISBI%sPhi
  alpha=0.0d0
  beta=0.0d0
  c=0.0d0
  DO i=2,nS
    DO j=1,i-1
      kFac=1
      DO k=1,pMax
        kFac=kfac*FLOAT(MAX(k-1,1))
        lFac=kFac
        DO l=1,sPhi(i)
          lFac=lFac*FLOAT(l+k-1)
          alpha(i,j,k)=alpha(i,j,k)+MISBI%a(i,j)%a(l)/lFac
        END DO
        kFac=kFac/MISBI%c(i)
      END DO
    END DO
  END DO
  DO i=1,ns
    DO j=1,i-1
      c(i)=c(i)+alpha(i,j,1)
    END DO  
  END DO  
  DO j=1,nS
    kFac=1
    DO k=1,pMax
      kFac=kfac*FLOAT(MAX(k-1,1))
      lFac=kFac
      DO l=1,sPhi(nS+1)
        lFac=lfac*FLOAT(l+k-1)
        beta(j,k)=beta(j,k)+MISBI%a(ns+1,j)%a(l)/lFac
      END DO
    END DO
  END DO

  Order=0.0d0
! Order 1
! N 
  DO i=1,ns
    Order(1)=Order(1)+beta(i,1)
  END DO  
  String(1)='.'
  Pos=FindTreeStringPos(String(1))
  pT=Pos(1)
  iT=Pos(2)
  WRITE(*,*) 'Order 1  1           N',Order(1),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0,TRIM(String(1))

! Order 2 
! N'N 
  DO i=1,ns
    Order(2)=Order(2)+beta(i,1)*c(i)
  END DO  
  String(2)='[.]'
  Pos=FindTreeStringPos(String(2))
  pT=Pos(1)
  iT=Pos(2)
  WRITE(*,*) "Order 2  1         N'N",Order(2),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),0.5d0,TRIM(String(2))
! LN 
  DO i=1,ns
    Order(3)=Order(3)+beta(i,2)
  END DO  
  String(3)='[o]'
  Pos=FindTreeStringPos(String(3))
  pT=Pos(1)
  iT=Pos(2)
  WRITE(*,*) 'Order 2  2          LN',Order(3),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),0.5d0,TRIM(String(3))

  IF (pMax>=3) THEN
!   Order 3 1  
!   N''(N,N)
    DO i=1,ns
      Order(4)=Order(4)+beta(i,1)*c(i)**2
    END DO  
    String(4)='[.,.]'
    Pos=FindTreeStringPos(String(4))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 3  1    N''(N,N)",Order(4),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/3.0d0,TRIM(String(4))
!   N'N'N
    DO i=1,ns
      DO j=1,i-1
        Order(5)=Order(5)+beta(i,1)*alpha(i,j,1)*c(j)
      END DO  
    END DO
    String(5)='[[.]]'
    Pos=FindTreeStringPos(String(5))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 3  2       N'N'N",Order(5),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/6.0d0,TRIM(String(5))
!   N'(LN)
    DO i=1,ns
      DO j=1,i-1
        Order(6)=Order(6)+beta(i,1)*alpha(i,j,2)
      END DO  
    END DO  
    String(6)='[[o]]'
    Pos=FindTreeStringPos(String(6))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 3  3      N'(LN)",Order(6),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/6.0d0,TRIM(String(6))
!   LN'N
    DO i=1,ns
      Order(7)=Order(7)+beta(i,2)*c(i)
    END DO  
    String(7)='[(.)]'
    Pos=FindTreeStringPos(String(7))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 3  4        LN'N",Order(7),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/6.0d0,TRIM(String(7))
!   LLN
    DO i=1,ns
      Order(8)=Order(8)+beta(i,3)
    END DO  
    String(8)='[(o)]'
    Pos=FindTreeStringPos(String(8))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 3  5         LLN",Order(8),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/6.0d0,TRIM(String(8))
  END IF
  IF (pMax>=4) THEN
!   N'''(N,N,N)   
    DO i=1,ns
      Order(9)=Order(9)+beta(i,1)*c(i)**3
    END DO  
    String(9)='[.,.,.]'
    Pos=FindTreeStringPos(String(9))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  1 N'''(N,N,N)",Order(9),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/4.0d0,TRIM(String(9))
!   N''(N',N,N)   
    DO i=1,ns
      DO j=1,i-1
        Order(10)=Order(10)+beta(i,1)*alpha(i,j,1)*c(i)*c(j)
      END DO  
    END DO  
    String(10)='[.,[.]]'
    Pos=FindTreeStringPos(String(10))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  2 N''(N',N,N)",Order(10),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/8.0d0,TRIM(String(10))
!   N''(L,N,N)   
    DO i=1,ns
      DO j=1,i-1
        Order(11)=Order(11)+beta(i,1)*alpha(i,j,2)*c(i)
      END DO  
    END DO  
    String(11)='[.,[o]]'
    Pos=FindTreeStringPos(String(11))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  3   N''(LN,N)",Order(11),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/8.0d0,TRIM(String(11))
!   N'N''(N,N)   
    DO i=1,ns
      DO j=1,i-1
        Order(12)=Order(12)+beta(i,1)*alpha(i,j,1)*c(j)**2
      END DO  
    END DO  
    String(12)='[[.,.]]'
    Pos=FindTreeStringPos(String(12))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  4  N'N''(N,N)",Order(12),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/12.0d0,TRIM(String(12))
!   LN''(N,N)   
    DO i=1,ns
      Order(13)=Order(13)+beta(i,2)*c(i)**2
    END DO  
    String(13)='[(.,.)]'
    Pos=FindTreeStringPos(String(13))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  5   LN''(N,N)",Order(13),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/12.0d0,TRIM(String(13))
!   N'N'N'N   
    DO i=1,ns
      DO j=1,i-1
        DO k=1,j-1
          Order(14)=Order(14)+beta(i,1)*alpha(i,j,1)*alpha(j,k,1)*c(k)
        END DO  
      END DO  
    END DO  
    String(14)='[[[.]]]'
    Pos=FindTreeStringPos(String(14))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  6     N'N'N'N",Order(14),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(14))
!   N'N'(LN)   
    DO i=1,ns
      DO j=1,i-1
        DO k=1,j-1
          Order(15)=Order(15)+beta(i,1)*alpha(i,j,1)*alpha(j,k,2)
        END DO  
      END DO  
    END DO  
    String(15)='[[[o]]]'
    Pos=FindTreeStringPos(String(15))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  7    N'N'(LN)",Order(15),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(15))
!   N'(LN'N)   
    DO i=1,ns
      DO j=1,i-1
        Order(16)=Order(16)+beta(i,1)*alpha(i,j,2)*c(j)
      END DO  
    END DO  
    String(16)='[[(.)]]'
    Pos=FindTreeStringPos(String(16))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  8    N'(LN'N)",Order(16),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(16))
!   N'(LLN)   
    DO i=1,ns
      DO j=1,i-1
        Order(17)=Order(17)+beta(i,1)*alpha(i,j,3)
      END DO  
    END DO  
    String(17)='[[(o)]]'
    Pos=FindTreeStringPos(String(17))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4  9     N'(LLN)",Order(17),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(17))
!   LN'N'N  
    DO i=1,ns
      DO j=1,i-1
        Order(18)=Order(18)+beta(i,2)*alpha(i,j,1)*c(j)
      END DO  
    END DO  
    String(18)='[([.])]'
    Pos=FindTreeStringPos(String(18))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4 10      LN'N'N",Order(18),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(18))
!   LN'(LN)  
    DO i=1,ns
      DO j=1,i-1
        Order(19)=Order(19)+beta(i,2)*alpha(i,j,2)
      END DO  
    END DO  
    String(19)='[([o])]'
    Pos=FindTreeStringPos(String(19))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4 11     LN'(LN)",Order(19),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(19))
!   LLN'N  
    DO i=1,ns
      Order(20)=Order(20)+beta(i,3)*c(i)
    END DO  
    String(20)='[((.))]'
    Pos=FindTreeStringPos(String(20))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4 12       LLN'N",Order(20),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(20))
!   LLLN  
    DO i=1,ns
      Order(21)=Order(21)+beta(i,4)
    END DO  
    String(21)='[((o))]'
    Pos=FindTreeStringPos(String(21))
    pT=Pos(1)
    iT=Pos(2)
    WRITE(*,*) "Order 4 12        LLLN",Order(21),Phi2MISBI(MISBI%nStage+1)%Phi(pT)%a(iT),1.0d0/24.0d0,TRIM(String(21))
  END IF

END SUBROUTINE OrderConditionsOwren
END MODULE CompareAlgorithm_Mod
