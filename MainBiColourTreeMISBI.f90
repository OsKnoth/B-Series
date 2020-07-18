PROGRAM TestTree
  USE MISBI_Mod
  USE BiColourTree_Mod
  USE BSeries_Mod
  USE CompareAlgorithm_Mod
  IMPLICIT NONE

  INTEGER :: p,pMaximum
  INTEGER :: i,i1,i2,j,k,l
  INTEGER :: NumberOrd
  INTEGER :: NumS,NumA
  REAL(8) :: Sym,Gam
  TYPE(TreeP_T), POINTER :: T1
  TYPE(TreeP_T), POINTER :: T2
  TYPE(STreeP_T), POINTER :: STree
  TYPE(MISBIMethod_T) :: MISBI
  TYPE(BSeries_T), ALLOCATABLE :: Phi2MISBI(:)
  TYPE(BPSeries_T), ALLOCATABLE :: Phi2PMISBI(:)
  INTEGER :: InputUnit=10
  CHARACTER(30) :: InputFileName
  CHARACTER(20) :: MISName='MISBIEB4'
  NAMELIST /Control/ pMaximum,    &
                     nStageSym,   &
                     dInsert,     &
                     gInsert,     &
                     nPhiSym,     &
                     MISName

  CALL get_command_argument(1,InputFileName)
  OPEN(UNIT=InputUnit,FILE=TRIM(InputFileName),STATUS='OLD')
  READ(InputUnit,NML=Control)
  CLOSE(InputUnit)

  CALL MISBIMethod(TRIM(MISName),MISBI,pMaximum)

  CALL InitTree(pMaximum)
  DO p=0,pMaximum
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      T1=>T1%Next
    END DO  
  END DO  

  CALL InitBSeries(pMaximum)

  DO p=0,pMaximum
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      CALL FindSTree(T1%TP)
      STree=>T1%TP%STree
      DO 
        IF (ASSOCIATED(STree)) THEN
          STree=>Stree%Next
        ELSE  
          EXIT
        END IF
      END DO
      T1=>T1%Next
    END DO  
  END DO  

  CALL MISBIBSeriesCompose(Phi2MISBI,Phi2PMISBI,MISBI,pMaximum)


  NumberOrd=1
  NumS=0
  NumA=0
  DO p=1,pMaximum
    WRITE(*,*) 'OrderCondition ',p
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      WRITE(*,*) '----------------------------------'
      WRITE(*,*) NumberOrd,T1%TP%Number,'Tree ',TRIM(T1%TP%TreeString)
      IF (INDEX(T1%TP%TreeString(1:1),'(')==0.AND. &
          INDEX(T1%TP%TreeString(1:1),'o')==0) THEN
        IF (INDEX(TRIM(T1%TP%TreeString),'(')==0.AND. &
            INDEX(TRIM(T1%TP%TreeString),'o')==0) THEN
          NumS=NumS+1  
          WRITE(*,*) NumS,'Classical Order ',p
          WRITE(*,*) TRIM(T1%TP%TreeString),E%Phi(p)%a(i),Phi2MISBI(MISBI%nStage+1)%Phi(p)%a(i) 
        ELSE  
          NumA=NumA+1  
          WRITE(*,*) NumA,'Additional Order ',p
          WRITE(*,*) TRIM(T1%TP%TreeString),E%Phi(p)%a(i),Phi2MISBI(MISBI%nStage+1)%Phi(p)%a(i) 
        END IF    
      ELSE  
        NumA=NumA+1  
        WRITE(*,*) NumA,'Additional Order ',p
        WRITE(*,*) TRIM(T1%TP%TreeString),E%Phi(p)%a(i),Phi2MISBI(MISBI%nStage+1)%Phi(p)%a(i) 
      END IF    
      IF (ABS(E%Phi(p)%a(i)-Phi2MISBI(MISBI%nStage+1)%Phi(p)%a(i))>1.d-8) THEN
        WRITE(*,*) 'OrderCondition failed'
      END IF  
      NumberOrd=NumberOrd+1
      T1=>T1%Next
    END DO  
  END DO  
  WRITE(*,*) 'OrderConditionsWenschKnoth '
  CALL OrderWenschKnoth(MISBI,Phi2MISBI,E,pMaximum)
  WRITE(*,*) 'OrderConditionsOwren '
  CALL OrderConditionsOwren(MISBI,Phi2MISBI,pMaximum)

END PROGRAM TestTree

