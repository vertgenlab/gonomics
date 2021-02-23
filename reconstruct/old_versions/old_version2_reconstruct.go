package reconstruct

/*
import (
// "fmt"
"github.com/vertgenlab/gonomics/fasta"
"github.com/vertgenlab/gonomics/dna"
)

type STree struct {
  	Name         string
  	BranchLength float64
    State         int
    Fasta       *fasta.Fasta
  	Left         *STree
  	Right        *STree

  }

  func Yhat(res []float64) int{
    var n float64
    var pos int
    for p,v:=range res {
        if v>n {

          n = v
          pos = p
        }

  }
  return pos
}

  func Get_branch(node *STree) []*STree {
    var branch []*STree
  if node.Left!=nil && node.Right !=nil{
  branch = append(branch, node)

  a:= Get_branch(node.Left)
  b:= Get_branch(node.Right)
  branch = append(branch, a...)
  branch= append(branch, b...)

  }

  return branch
  }


  func Get_leaf(node *STree) []*STree{
    var leaf [] *STree
    if node.Left !=nil && node.Right!=nil{
    a:= Get_leaf(node.Left)
    b:= Get_leaf(node.Right)
    leaf = append(leaf, a...)
    leaf= append(leaf, b...)}

    if node.Left==nil && node.Right==nil{
      leaf = append(leaf, node)

    }
    return leaf
  }

  func Get_branch_fixed(root *STree, fixed *STree) []*STree {
    var node *STree
    branches:= Get_branch(root)
    var branches_new []*STree

  for num:=0; num<len(branches); num++{
  node =branches[num]
  if node != fixed{
  branches_new = append(branches_new, node)

  }}


  return branches_new
  }


func Tp(node *STree) float64 {
  prod:= 1.0
if node.Left != nil{
  a:= node.State
  b:= node.Left.State
  t:= node.Left.BranchLength

  var p float64
  switch{
    case a > 3 || b >3:
      p = 0
    case a == b:
      p = 1-t
    default:
      p = t * .25
}
prod = prod * p
l:= Tp(node.Left)
prod = prod * l
}


if node.Right != nil{
  a:= node.State
  b:= node.Right.State
  t:= node.Right.BranchLength

  var p float64
  switch{
    case a > 3 || b >3:
      p = 0
    case a == b:
      p = 1-t
    default:
      p = t * .25
}
prod = prod * p
r:= Tp(node.Right)
prod = prod *r
}

return prod
}

func Numerator(root *STree, fixed *STree, num int) float64 {
  sum:= 0.0
var node *STree
var ans float64
  branches:= Get_branch_fixed(root, fixed)
node = branches[num]

  for i:= 0; i<4; i++{
node.State = i


if len(branches)-num >=0 && len(branches)-num < len(branches){
a:= Numerator(root,fixed, num-1)
sum+=a}

if node.Name == branches[0].Name{
ans = Tp(root)
sum += ans

}



  }



  return sum

  }

  func Denominator(root *STree, num int) float64 {
  sum:= 0.0

  var ans float64

  branches:= Get_branch(root)



  if len(branches)-num >=0 && len(branches)-num < len(branches){
  for i:= 0; i<4; i++{

  branches[len(branches)-num].State = i


  if branches[len(branches)-num].Name == branches[len(branches)-1].Name{
  ans = Tp(root)
  sum += ans

  }

  a:= Denominator(root, num-1)
  sum+=a
  }

  }
  return sum

  }

func Eqn (root *STree) []int  {
var states [] int
internal_nodes := Get_branch(root)
  denominator:= .25 *Denominator(root, len(internal_nodes))
var node *STree
for i:= 0; i<len(internal_nodes); i++{
  node = internal_nodes[i]
  var res [] float64
  var yhat int
for j:= 0; j<4; j++{
  node.State = j
  sum:= .25 *Numerator(root, node, len(internal_nodes)-2)
  res = append(res, sum/denominator)

}
yhat = Yhat(res)
// fmt.Print(node.Name, " = ", yhat, " with probability ",res[yhat] , "\n")
states = append(states, yhat)
}
return states
}




func Loop_sites (node *STree) []fasta.Fasta{
leafs:=Get_leaf(node)
branches := Get_branch(node)
var fastas []*fasta.Fasta
for j:=0; j<len(leafs[0].Fasta.Seq); j++{
  for i :=0; i<len(leafs); i++{

    leafs[i].State = int(leafs[i].Fasta.Seq[j])
  }

yhat:= Eqn(node)
for i := 0; i <len(branches); i++{

branches[i].Fasta.Seq = append(branches[i].Fasta.Seq, []dna.Base{dna.Base(yhat[i])}...)

}}

for i :=0; i <len(branches); i++{

  fastas=append(fastas, branches[i].Fasta)
}

var fas []fasta.Fasta
for i:=0; i<len(fastas); i++{
fas = append(fas, *fastas[i])

}

return fas
}
*/
