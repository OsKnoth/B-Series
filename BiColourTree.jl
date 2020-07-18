module BiColourTree

#TYPE Tree_T
#    INTEGER :: Order=0
#    INTEGER :: Number=0
#    INTEGER :: NumberOrder=0
#    TYPE(TreeP_T), POINTER :: Childs(:)=>NULL()
#    INTEGER :: NumberChilds=0
#    INTEGER, POINTER :: OrderChilds(:)=>NULL()
#    CHARACTER(LenTreeString) :: TreeString=''
#    CHARACTER :: Colour
#    LOGICAL :: NoWhiteLeaf=.TRUE.
#    INTEGER :: LenSTree=0
#    TYPE(STreeP_T), POINTER :: STree=>NULL()
#    TYPE(TreeP_T), POINTER :: Powers(:)=>NULL()
#    INTEGER :: NumberPowers=0
#  END TYPE Tree_T

mutable struct Tree_T
  Order::Int 
  Number::Int 
  NumberOrder::Int
  Child::Array{Tree_T,1}
  NumberChilds::Int
  OrderChilds::Array{Int,1}
end


function Tree_T(Order::Int,Number::Int)
  NumberOrder=0
  Child=Array{Tree_T,1}()
  Tree_T(Order,Number,NumberOrder,Child)
end  
function Tree_T(Order::Int,Number::Int,Child::Array{Tree_T,1})
  NumberOrder=0
  Tree_T(Order,Number,NumberOrder,Child)
end  
#if isless(Order, 0) && throw(DomainError("Order must be nonnegative"))


#Tree1=Tree_T(21,[]::Array{Int,1})
Tree1=Tree_T(2,1)
Tree2=Tree_T(2,2)

Tree3=Tree_T(1,2,[Tree1,Tree2])

#println("Order ", Tree1.Order)
#println("Order ", Tree1.NumberOrder)

end
