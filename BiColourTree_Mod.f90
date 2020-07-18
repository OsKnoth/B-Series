MODULE BiColourTree_Mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: LenTreeString=200
  TYPE Tree_T
    INTEGER :: Order=0
    INTEGER :: Number=0
    INTEGER :: NumberOrder=0
    TYPE(TreeP_T), POINTER :: Childs(:)=>NULL()
    INTEGER :: NumberChilds=0
    INTEGER, POINTER :: OrderChilds(:)=>NULL()
    CHARACTER(LenTreeString) :: TreeString=''
    CHARACTER :: Colour
    LOGICAL :: NoWhiteLeaf=.TRUE.
    INTEGER :: LenSTree=0
    TYPE(STreeP_T), POINTER :: STree=>NULL()
!   INTEGER, POINTER :: OrderSTree(:)=>NULL()
    TYPE(TreeP_T), POINTER :: Powers(:)=>NULL()
    INTEGER :: NumberPowers=0
  END TYPE Tree_T

  TYPE Forest_T
    INTEGER :: LenForest=0
    TYPE(TreeP_T), POINTER :: Trees(:)=>NULL()
    INTEGER, POINTER :: OrderTrees(:)=>NULL()
  END TYPE Forest_T

  TYPE TreeP_T
    TYPE(Tree_T), POINTER :: TP=>NULL()
    TYPE(TreeP_T), POINTER :: Next=>NULL()
  END TYPE TreeP_T

  TYPE STreeP_T
    TYPE(Tree_T), POINTER :: TP=>NULL()
    INTEGER, POINTER :: ListCuts(:)
    INTEGER, POINTER :: NumberListCuts(:)
    TYPE(STreeP_T), POINTER :: Next=>NULL()
  END TYPE STreeP_T

  TYPE ListTree_T 
    INTEGER :: LenListTree=0
    TYPE(TreeP_T), POINTER :: T=>NULL()
  END TYPE ListTree_T  

  TYPE(ListTree_T), ALLOCATABLE, TARGET :: ListTree(:)
  INTEGER, PRIVATE :: pMax
  TYPE(TreeP_T), ALLOCATABLE :: NumberToTree(:)

  CHARACTER(LenTreeString) :: TreeString
  INTEGER :: PosTreeString
  INTEGER :: NumberTrees


CONTAINS

FUNCTION FindTreeStringPos(TreeString)
  INTEGER :: FindTreeStringPos(2)
  CHARACTER(LenTreeString) :: TreeString
  INTEGER :: i,p
  TYPE(TreeP_T), POINTER :: T1
  S1:DO p=0,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      IF (T1%TP%TreeString==TreeString) THEN
        FindTreeStringPos(1)=p
        FindTreeStringPos(2)=i
        EXIT S1
      END IF  
      T1=>T1%Next
    END DO  
  END DO S1  
END FUNCTION FindTreeStringPos

SUBROUTINE InitTree(pMaximum)
  INTEGER :: pMaximum

  INTEGER :: i,i1,i2,p,TreeNumber
  TYPE(TreeP_T), POINTER :: T
  TYPE(TreeP_T), POINTER :: T1
  TYPE(TreeP_T), POINTER :: T2
  TYPE(TreeP_T) :: TWhite

  pMax=pMaximum
  TreeNumber=0
  ALLOCATE(ListTree(0:pMax))
  DO i=0,pMax
    ListTree(i)%LenListTree=0
  END DO  

  ALLOCATE(ListTree(0)%T)
  ALLOCATE(ListTree(0)%T%TP)
  ListTree(0)%T%TP%Order=0
  TreeNumber=TreeNumber+1
  ListTree(0)%T%TP%Number=TreeNumber

  ListTree(0)%T%TP%NumberOrder=1
  ListTree(0)%LenListTree=1
  CALL TreeToString(ListTree(0)%T%TP)

  ListTree(1)%LenListTree=2
  ALLOCATE(ListTree(1)%T)
  T=>ListTree(1)%T
  ALLOCATE(T%TP)
  T%TP%Order=1
  TreeNumber=TreeNumber+1
  T%TP%Number=TreeNumber
  T%TP%NumberOrder=1
  T%TP%Colour='b'
  T%TP%NoWhiteLeaf=.TRUE.
  CALL TreeToString(T%TP)
  ALLOCATE(T%Next)
  T=>T%Next
  ALLOCATE(T%TP)
  T%TP%Order=1
  TreeNumber=TreeNumber+1
  T%TP%Number=TreeNumber
  T%TP%NumberOrder=2
  T%TP%Colour='w'
  T%TP%NoWhiteLeaf=.FALSE.
  CALL TreeToString(T%TP)

  ALLOCATE(TWhite%TP)
  TWhite%TP%Order=1
  TWhite%TP%Colour='w'
  CALL TreeToString(TWhite%TP)

  DO p=1,pMax-1
    DO i=1,p
      T1=>ListTree(i)%T
      DO i1=1,ListTree(i)%LenListTree
        T2=>ListTree(p+1-i)%T
        DO i2=1,ListTree(p+1-i)%LenListTree
          CALL Prod(T1%TP,T2%TP)
          T2=>T2%Next
        END DO  
        T1=>T1%Next
      END DO  
    END DO  
!   T2=>ListTree(p)%T
!   DO i1=1,ListTree(p)%LenListTree
!     CALL Prod(TWhite%TP,T2%TP)
!     T2=>T2%Next
!   END DO  
    T1=>ListTree(p+1)%T
    DO i=1,ListTree(p+1)%LenListTree
      TreeNumber=TreeNumber+1
      T1%TP%Number=TreeNumber
      T1%TP%NumberOrder=i
      T1=>T1%Next
    END DO  
  END DO  
  NumberTrees=TreeNumber
  ALLOCATE(NumberToTree(TreeNumber))
  TreeNumber=1
  DO p=0,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      NumberToTree(TreeNumber)%TP=>T1%TP
      TreeNumber=TreeNumber+1
      T1=>T1%Next
    END DO  
  END DO  
  
END SUBROUTINE InitTree

SUBROUTINE PowersTree(T)
  TYPE(Tree_T), POINTER :: T

  INTEGER :: i,Order
  TYPE(Tree_T) :: P
  CHARACTER(LenTreeString) :: TreeString

  Order=T%Order
  T%NumberPowers=pMax-Order
  ALLOCATE(T%Powers(0:T%NumberPowers))
  TreeString=T%TreeString
  T%Powers(0)%TP=>NumberToTree(FindTreeString(TreeString))%TP
  DO i=1,T%NumberPowers
    IF (TreeString=='') THEN
      TreeString='o'
    ELSE  
      TreeString='('//TRIM(TreeString)//')'
    END IF  
    T%Powers(i)%TP=>NumberToTree(FindTreeString(TreeString))%TP
  END DO
END SUBROUTINE PowersTree


FUNCTION CompareTree(T1,T2)
  INTEGER :: CompareTree
  TYPE(Tree_T), POINTER :: T1
  TYPE(Tree_T), POINTER :: T2
  
  IF (T1%Number<T2%Number) THEN
    CompareTree=-1
  ELSE IF (T1%Number==T2%Number) THEN   
    CompareTree=0
  ELSE   
     CompareTree=1
  END IF   
END FUNCTION CompareTree

SUBROUTINE InsertTree(T,Forest,Insert)
  TYPE(Tree_T), POINTER :: T
  TYPE(TreeP_T), POINTER :: Forest
  LOGICAL :: Insert

  TYPE(TreeP_T), POINTER :: Current
  TYPE(TreeP_T), POINTER :: Previous
  TYPE(TreeP_T), POINTER :: TreeP

  Current=>Forest
  IF (CompareTree(T,Current%TP)<0) THEN
    ALLOCATE(TreeP)
    TreeP%TP=>T
    Forest=>TreeP
    Forest%Next=>Current
    Insert=.TRUE.
  ELSE  
    DO
      IF (ASSOCIATED(Current)) THEN
        IF (CompareTree(T,Current%TP)<0) THEN
          ALLOCATE(TreeP)
          TreeP%TP=>T
          Previous%Next=>TreeP
          TreeP%Next=>Current
          Insert=.TRUE.
          EXIT
        ELSE IF (CompareTree(T,Current%TP)==0) THEN
          Insert=.FALSE.
          EXIT
        END IF  
      ELSE  
        ALLOCATE(TreeP)
        TreeP%TP=>T
        Current%Next=>TreeP
        Insert=.TRUE.
        EXIT
      END IF
      Previous=>Current
      Current=>Current%Next
    END DO  
  END IF  
END SUBROUTINE InsertTree

FUNCTION ContainTree(T,Forest)
  LOGICAL :: ContainTree
  TYPE(Tree_T), POINTER :: T
  TYPE(TreeP_T), POINTER :: Forest

  TYPE(TreeP_T), POINTER :: Current

  Current=>Forest
  DO
  END DO
END FUNCTION ContainTree

SUBROUTINE SortTree(T)
  TYPE(Tree_T), POINTER :: T

  INTEGER :: i,j,iTemp
  INTEGER :: NewNumberChilds
  LOGICAL :: InsertChild(T%NumberChilds)
  TYPE(Tree_T), POINTER :: TreeTemp
  TYPE(TreeP_T), POINTER :: Childs(:)
  INTEGER, POINTER :: OrderChilds(:)

  NewNumberChilds=T%NumberChilds 
  InsertChild=.TRUE.
  DO i=1,T%NumberChilds
    IF (InsertChild(i)) THEN
      DO j=i+1,T%NumberChilds
        IF (T%Childs(i)%TP%Number==T%Childs(j)%TP%Number) THEN
          T%OrderChilds(i)=T%OrderChilds(i)+1
          NewNumberChilds=NewNumberChilds-1
          InsertChild(j)=.FALSE.
        ELSE IF (T%Childs(i)%TP%Number>T%Childs(j)%TP%Number.AND.InsertChild(j)) THEN  
          TreeTemp=>T%Childs(i)%TP
          T%Childs(i)%TP=>T%Childs(j)%TP
          T%Childs(j)%TP=>TreeTemp
          iTemp=T%OrderChilds(i)
          T%OrderChilds(i)=T%OrderChilds(j)
          T%OrderChilds(j)=iTemp
        END IF  
      END DO  
    END IF  
  END DO  
  IF (NewNumberChilds<T%NumberChilds) THEN
    ALLOCATE(Childs(NewNumberChilds))
    ALLOCATE(OrderChilds(NewNumberChilds))
    iTemP=1
    DO i=1,T%NumberChilds
      IF (InsertChild(i)) THEN
        Childs(iTemp)%TP=>T%Childs(i)%TP
        OrderChilds(iTemp)=T%OrderChilds(i)
        iTemp=iTemp+1
      END IF  
    END DO  
    DEALLOCATE(T%Childs)
    DEALLOCATE(T%OrderChilds)
    T%Childs=>Childs
    T%OrderChilds=>OrderChilds
    T%NumberChilds=NewNumberChilds
  END IF  
  CALL TreeToString(T)
  CALL FindNumber(T)
END SUBROUTINE SortTree

SUBROUTINE FindNumber(T)
  TYPE(Tree_T), POINTER :: T

  INTEGER :: i,p
  TYPE(TreeP_T), POINTER :: T1

  S1:DO p=0,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      IF (T1%TP%TreeString==T%TreeString) THEN
        T%Number=T1%TP%Number
        EXIT S1
      END IF  
      T1=>T1%Next
    END DO  
  END DO S1  
END SUBROUTINE FindNumber

FUNCTION FindTreeString(TreeString)
  INTEGER :: FindTreeString
  CHARACTER(LenTreeString) :: TreeString
  INTEGER :: i,p
  TYPE(TreeP_T), POINTER :: T1
  S1:DO p=0,pMax
    T1=>ListTree(p)%T
    DO i=1,ListTree(p)%LenListTree
      IF (T1%TP%TreeString==TreeString) THEN
        FindTreeString=T1%TP%Number
        EXIT S1
      END IF  
      T1=>T1%Next
    END DO  
  END DO S1  
END FUNCTION FindTreeString

  
SUBROUTINE FindSTree(T)
  TYPE(Tree_T), POINTER :: T

  INTEGER :: i,iC,j,k,l,l1,l2,l1End,l2End,m,LenSTree,LenSTreeRed
  INTEGER :: NumberChilds
  INTEGER :: Order
  INTEGER :: LenListCuts
  INTEGER :: TreeList(NumberTrees)
  INTEGER, POINTER :: NumberList(:)
  INTEGER, POINTER :: NumberListOrder(:)
  INTEGER, POINTER :: ListCut(:,:)
  TYPE(STreeP_T), POINTER :: STree(:)
  TYPE(STreeP_T), POINTER :: STreeTemp
  TYPE(STreeP_T), POINTER :: Current

  IF (T%Order==0) THEN
    T%LenSTree=1
    ALLOCATE(T%STree)
    T%STree%TP=>T
    ALLOCATE(T%STree%ListCuts(1))
    ALLOCATE(T%STree%NumberListCuts(1))
    T%STree%ListCuts(1)=1
    T%STree%NumberListCuts(1)=1
  ELSE IF (T%Order==1) THEN
    T%LenSTree=2
    ALLOCATE(T%STree)
    T%STree%TP=>NumberToTree(1)%TP
    ALLOCATE(T%STree%ListCuts(1))
    ALLOCATE(T%STree%NumberListCuts(1))
    T%STree%ListCuts(1)=T%Number
    T%STree%NumberListCuts(1)=1
    ALLOCATE(T%STree%Next)
    Current=>T%STree%Next
    Current%TP=>T
    ALLOCATE(Current%ListCuts(1))
    ALLOCATE(Current%NumberListCuts(1))
    Current%ListCuts(1)=1
    Current%NumberListCuts(1)=1
  ELSE  
    LenSTree=1
    NumberChilds=0
    DO i=1,T%NumberChilds
      DO j=1,T%OrderChilds(i)
        NumberChilds=NumberChilds+1
        LenSTree=LenSTree*T%Childs(i)%TP%LenSTree
      END DO  
    END DO  
    ALLOCATE(STree(LenSTree))
    ALLOCATE(ListCut(LenSTree,NumberTrees))
    ListCut=0
    DO k=1,LenSTree
      ALLOCATE(STree(k)%TP)
      ALLOCATE(STree(k)%TP%Childs(NumberChilds))
      ALLOCATE(STree(k)%TP%OrderChilds(NumberChilds))
      STree(k)%TP%NumberChilds=NumberChilds
      STree(k)%TP%OrderChilds=1
      STree(k)%TP%Colour=T%Colour
    END DO  
    l1End=LenSTree
    l2End=1
    iC=0
    DO i=1,T%NumberChilds
      DO l=1,T%OrderChilds(i)
        iC=iC+1
        l1End=l1End/T%Childs(i)%TP%LenSTree
        k=1
        DO l1=1,l1End 
          STreeTemp=>T%Childs(i)%TP%STree
          DO j=1,T%Childs(i)%TP%LenSTree
            DO l2=1,l2end 
              STree(k)%TP%Childs(iC)%TP=>STreeTemp%TP
!             ListCut(k,STreeTemp%TP%Number)=1
              DO m=1,SIZE(STreeTemp%ListCuts)
                ListCut(k,STreeTemp%ListCuts(m))=ListCut(k,STreeTemp%ListCuts(m))+STreeTemp%NumberListCuts(m)
              END DO  
              k=k+1
            END DO  
            STreeTemp=>STreeTemp%Next 
          END DO  
        END DO  
        l2End=l2End*T%Childs(i)%TP%LenSTree
      END DO  
    END DO  
    ALLOCATE(NumberList(LenSTree))

    LenSTreeRed=0
    DO k=1,LenSTree
      Order=1
      DO i=1,STree(k)%TP%NumberChilds
        Order=Order+STree(k)%TP%Childs(i)%TP%Order
      END DO  
      STree(k)%TP%Order=Order
      CALL SortTree(STree(k)%TP)
      NumberList(k)=STree(k)%TP%Number
      IF (STree(k)%TP%Number>0) THEN
        LenSTreeRed=LenSTreeRed+1
      END IF  
      DEALLOCATE(STree(k)%TP%Childs)
      DEALLOCATE(STree(k)%TP%OrderChilds)
      DEALLOCATE(STree(k)%TP)
    END DO  
!   CALL CompressList(NumberListRed)
    DEALLOCATE(STree)
    ALLOCATE(T%STree)
    T%STree%TP=>NumberToTree(1)%TP
    ALLOCATE(T%STree%ListCuts(1))
    T%STree%ListCuts(1)=T%Number
    ALLOCATE(T%STree%NumberListCuts(1))
    T%STree%NumberListCuts(1)=1
    Current=>T%STree
!   ALLOCATE(Current%Next)
!   Current=>Current%Next
!   Current%TP=>NumberToTree(T%Number)%TP
!   ALLOCATE(Current%ListCuts(1))
!   Current%ListCuts(1)=1
!   ALLOCATE(Current%NumberListCuts(1))
!   Current%NumberListCuts(1)=1
    LenSTree=1
    DO k=1,SIZE(NumberList)
      IF (NumberList(k)>0) THEN
        LenSTree=LenSTree+1
        ALLOCATE(Current%Next)
        Current=>Current%Next
        Current%TP=>NumberToTree(NumberList(k))%TP
        LenListCuts=0
        DO j=2,NumberTrees
          IF (ListCut(k,j)>0) THEN
            LenListCuts=LenListCuts+1
          END IF  
        END DO  
        IF (LenListCuts==0) THEN
          LenListCuts=1
        END IF  
        ALLOCATE(Current%ListCuts(LenListCuts))
        ALLOCATE(Current%NumberListCuts(LenListCuts))
        LenListCuts=0
        DO j=2,NumberTrees
          IF (ListCut(k,j)>0) THEN
            LenListCuts=LenListCuts+1
            Current%ListCuts(LenListCuts)=j
            Current%NumberListCuts(LenListCuts)=ListCut(k,j)
          END IF  
        END DO  
        IF (LenListCuts==0) THEN
          Current%ListCuts(1)=1
          Current%NumberListCuts(1)=1
        END IF  
      END IF
    END DO  
    T%LenSTree=LenSTree
    ALLOCATE(NumberListOrder(MAXVAL(NumberList)))
    NumberListOrder=0
    DO k=1,SIZE(NumberList)
      IF (NumberList(k)>0) THEN
        NumberListOrder(NumberList(k))=NumberListOrder(NumberList(k))+1
      END IF
    END DO
    DEALLOCATE(NumberList)
    DEALLOCATE(NumberListOrder)
    DEALLOCATE(ListCut)
  END IF  

END SUBROUTINE FindSTree

RECURSIVE SUBROUTINE PrintTree(T,Ref)
  TYPE(Tree_T), POINTER :: T
  INTEGER, OPTIONAL :: Ref

  INTEGER :: i,j,RefLoc

  IF (.NOT.PRESENT(Ref)) THEN
    RefLoc=1
    PosTreeString=1
  ELSE
    RefLoc=Ref
  END IF  
  IF (T%Order>0) THEN
    IF (T%NumberChilds>0) THEN
      TreeString(PosTreeString:PosTreeString)='['
      PosTreeString=PosTreeString+1
      DO i=1,T%NumberChilds
        DO j=1,T%OrderChilds(i)
          CALL PrintTree(T%Childs(i)%TP,RefLoc+1)
          TreeString(PosTreeString:PosTreeString)=','
          PosTreeString=PosTreeString+1
        END DO
        PosTreeString=PosTreeString-1
      END DO  
      TreeString(PosTreeString:PosTreeString)=']'
      PosTreeString=PosTreeString+1
    ELSE
      TreeString(PosTreeString:PosTreeString)='o'
      PosTreeString=PosTreeString+1
    END IF   
  END IF
  IF (RefLoc==1) THEN
    WRITE(*,*)  TreeString(1:PosTreeString-1)
  END IF  
END SUBROUTINE PrintTree 

RECURSIVE SUBROUTINE TreeToString(T,Ref)
  TYPE(Tree_T), POINTER :: T
  INTEGER, OPTIONAL :: Ref

  INTEGER :: i,j,RefLoc

  IF (.NOT.PRESENT(Ref)) THEN
    RefLoc=1
    PosTreeString=1
  ELSE
    RefLoc=Ref
  END IF  
  IF (T%Order>0) THEN
    IF (T%NumberChilds>1) THEN
      IF (T%Colour=='w') THEN
        TreeString(PosTreeString:PosTreeString)='('
      ELSE  
        TreeString(PosTreeString:PosTreeString)='['
      END IF  
      PosTreeString=PosTreeString+1
      DO i=1,T%NumberChilds
        IF (T%Childs(i)%TP%Order>0) THEN
          DO j=1,T%OrderChilds(i)
            CALL TreeToString(T%Childs(i)%TP,RefLoc+1)
            TreeString(PosTreeString:PosTreeString)=','
            PosTreeString=PosTreeString+1
          END DO
          PosTreeString=PosTreeString-1
        END IF
        IF (i<T%NumberChilds) THEN
          TreeString(PosTreeString:PosTreeString)=','
          PosTreeString=PosTreeString+1
        END IF
      END DO  
      IF (T%Colour=='w') THEN
        TreeString(PosTreeString:PosTreeString)=')'
      ELSE  
        TreeString(PosTreeString:PosTreeString)=']'
      END IF  
      PosTreeString=PosTreeString+1
    ELSE IF (T%NumberChilds==1) THEN
      IF (T%Childs(1)%TP%Number>1) THEN
        IF (T%Colour=='w') THEN
          TreeString(PosTreeString:PosTreeString)='('
        ELSE  
          TreeString(PosTreeString:PosTreeString)='['
        END IF  
        PosTreeString=PosTreeString+1
        DO i=1,T%NumberChilds
          DO j=1,T%OrderChilds(i)
            CALL TreeToString(T%Childs(i)%TP,RefLoc+1)
            TreeString(PosTreeString:PosTreeString)=','
            PosTreeString=PosTreeString+1
          END DO
          PosTreeString=PosTreeString-1
        END DO  
        IF (T%Colour=='w') THEN
          TreeString(PosTreeString:PosTreeString)=')'
        ELSE  
          TreeString(PosTreeString:PosTreeString)=']'
        END IF  
        PosTreeString=PosTreeString+1
      ELSE  
        IF (T%Order>0) THEN
          IF (T%Colour=='w') THEN
            TreeString(PosTreeString:PosTreeString)='o'
          ELSE  
            TreeString(PosTreeString:PosTreeString)='.'
          END IF  
          PosTreeString=PosTreeString+1
        END IF  
      END IF  
    ELSE
      IF (T%Colour=='w') THEN
        TreeString(PosTreeString:PosTreeString)='o'
      ELSE  
        TreeString(PosTreeString:PosTreeString)='.'
      END IF  
      PosTreeString=PosTreeString+1
    END IF   
  END IF   
  IF (RefLoc==1) THEN
    T%TreeString(1:PosTreeString-1)=TreeString(1:PosTreeString-1)
  END IF  
END SUBROUTINE TreeToString 

RECURSIVE FUNCTION Gamma(T,Ref) RESULT(Gam)
  REAL(8) :: Gam
  TYPE(Tree_T), POINTER :: T
  INTEGER, OPTIONAL :: Ref

  INTEGER :: i,j,RefLoc

  Gam=MAX(FLOAT(T%Order),1.0d0)
  IF (.NOT.PRESENT(Ref)) THEN
    RefLoc=1
  ELSE
    RefLoc=Ref
  END IF  
  IF (T%NumberChilds>0.AND.T%Order>0) THEN
    DO i=1,T%NumberChilds
      Gam=Gam*Gamma(T%Childs(i)%TP,RefLoc+1)**T%OrderChilds(i)
    END DO  
  END IF   
END FUNCTION Gamma 

RECURSIVE FUNCTION Symmetry(T,Ref) RESULT(Sym)
  REAL(8) :: Sym
  TYPE(Tree_T), POINTER :: T
  INTEGER, OPTIONAL :: Ref

  INTEGER :: i,j,RefLoc

  Sym=1.0d0
  IF (.NOT.PRESENT(Ref)) THEN
    RefLoc=1
  ELSE
    RefLoc=Ref
  END IF  
  IF (T%NumberChilds>0.AND.T%Order>0) THEN
    DO i=1,T%NumberChilds
      Sym=Sym*Fak(T%OrderChilds(i)) &
         *Symmetry(T%Childs(i)%TP,RefLoc+1)**T%OrderChilds(i)
    END DO  
  END IF   
END FUNCTION Symmetry 

FUNCTION Fak(n)
  INTEGER :: Fak
  INTEGER :: n

  INTEGER :: i

  Fak=1
  DO i=1,n
    Fak=Fak*i
  END DO  
END FUNCTION Fak

SUBROUTINE Prod(T1,T2)
  TYPE(Tree_T), POINTER :: T1,T2

  INTEGER :: i,j,jj
  INTEGER :: T2Order
  LOGICAL :: Append
  TYPE(TreeP_T), POINTER :: Current,Previous
  TYPE(TreeP_T), POINTER :: TProd=>NULL()

! Check iff T2 is a child of T1

  Append=.TRUE.
  ALLOCATE(TProd)
  ALLOCATE(TProd%TP)
  TProd%TP%Order=T1%Order+T2%Order
  TProd%TP%Colour=T1%Colour
  TProd%TP%NoWhiteLeaf=T1%NoWhiteLeaf
! IF (.NOT.T2%NoWhiteLeaf) THEN
!   TProd%TP%NoWhiteLeaf=.FALSE.
! END IF  
  DO i=1,T1%NumberChilds
    IF (T2%Order==T1%Childs(i)%TP%Order.AND. &
        T2%Number==T1%Childs(i)%TP%Number) THEN
      TProd%TP%NumberChilds=T1%NumberChilds
      ALLOCATE(TProd%TP%Childs(TProd%TP%NumberChilds))
      DO j=1,TProd%TP%NumberChilds
        TProd%TP%Childs(j)%TP=>T1%Childs(j)%TP
      END DO  
      ALLOCATE(TProd%TP%OrderChilds(TProd%TP%NumberChilds))
      TProd%TP%OrderChilds=T1%OrderChilds
      TProd%TP%OrderChilds(i)=TProd%TP%OrderChilds(i)+1
      Append=.FALSE.
      EXIT
    END IF  
  END DO  
  IF (Append) THEN
    TProd%TP%NumberChilds=T1%NumberChilds+1
    ALLOCATE(TProd%TP%Childs(TProd%TP%NumberChilds))
    ALLOCATE(TProd%TP%OrderChilds(TProd%TP%NumberChilds))
    TProd%TP%OrderChilds(TProd%TP%NumberChilds)=0
    jj=1
    IF (TProd%TP%NumberChilds>1) THEN
      T2Order=T2%Order
      DO j=1,TProd%TP%NumberChilds-1
        IF (T1%Childs(j)%TP%Order<T2Order) THEN
          TProd%TP%Childs(jj)%TP=>T1%Childs(j)%TP
          TProd%TP%OrderChilds(jj)=T1%OrderChilds(j)
          jj=jj+1
        ELSE IF (T1%Childs(j)%TP%Order==T2Order) THEN  
          IF (T1%Childs(j)%TP%Number<T2%Number) THEN
            TProd%TP%Childs(jj)%TP=>T1%Childs(j)%TP
            TProd%TP%OrderChilds(jj)=T1%OrderChilds(j)
            jj=jj+1
          ELSE  
            TProd%TP%Childs(jj)%TP=>T2
            TProd%TP%OrderChilds(jj)=1
            T2Order=TProd%TP%Order
            jj=jj+1
            TProd%TP%Childs(jj)%TP=>T1%Childs(j)%TP
            TProd%TP%OrderChilds(jj)=T1%OrderChilds(j)
            jj=jj+1
          END IF  
        ELSE
          TProd%TP%Childs(jj)%TP=>T2
          TProd%TP%OrderChilds(jj)=1
          T2Order=TProd%TP%Order
          jj=jj+1
          TProd%TP%Childs(jj)%TP=>T1%Childs(j)%TP
          TProd%TP%OrderChilds(jj)=T1%OrderChilds(j)
          jj=jj+1
        END IF  
      END DO  
      IF (TProd%TP%OrderChilds(TProd%TP%NumberChilds)==0) THEN
        TProd%TP%Childs(TProd%TP%NumberChilds)%TP=>T2
        TProd%TP%OrderChilds(TProd%TP%NumberChilds)=1
      END IF  
    ELSE
      TProd%TP%Childs(1)%TP=>T2
      TProd%TP%OrderChilds(1)=1
    END IF
  END IF
  CALL TreeToString(TProd%TP)
  IF (ASSOCIATED(ListTree(TProd%TP%Order)%T)) THEN
    jj=0
    Current=>ListTree(TProd%TP%Order)%T
    DO 
      jj=jj+1
      IF (TProd%TP%TreeString<Current%TP%TreeString) THEN
        TProd%Next=>Current
        IF (jj==1) THEN
          ListTree(TProd%TP%Order)%T=>TProd
        ELSE
          Previous%Next=>TProd
        END IF  
        ListTree(TProd%TP%Order)%LenListTree=ListTree(TProd%TP%Order)%LenListTree+1
        EXIT
      ELSE IF (TProd%TP%TreeString==Current%TP%TreeString) THEN
        DEALLOCATE(TProd%TP%Childs)
        DEALLOCATE(TProd%TP%OrderChilds)
        DEALLOCATE(TProd)
        EXIT
      END IF  
      IF (ASSOCIATED(Current%Next)) THEN
        Previous=>Current
        Current=>Current%Next
      ELSE  
        Current%Next=>TProd
        ListTree(TProd%TP%Order)%LenListTree=ListTree(TProd%TP%Order)%LenListTree+1
        EXIT
      END IF  
    END DO
  ELSE  
    ListTree(TProd%TP%Order)%T=>TProd
    ListTree(TProd%TP%Order)%LenListTree=1
  END IF  
END SUBROUTINE Prod

SUBROUTINE CompressList(List)

  INTEGER, POINTER :: List(:)
  INTEGER :: i,j,iList,Member
  INTEGER :: TempList(SIZE(List))
  LOGICAL :: Insert

  iList=0
  S1:DO i=1,SIZE(List)
    Member=List(i) 
    Insert=.TRUE.
    S2:DO j=1,iList
      IF (Member==TempList(j)) THEN
        Insert=.FALSE.
        EXIT S2
      END IF   
    END DO S2  
    IF (Insert) THEN
      iList=iList+1
      TempList(iList)=Member
    END IF  
  END DO S1
  DEALLOCATE(List)
  ALLOCATE(List(1:iList))
  List=TempList(1:iList)

END SUBROUTINE CompressList

END MODULE BiColourTree_Mod
