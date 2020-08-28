package expandedTree

import (
	"bufio"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"os"
	"strings"
)

//Tree structure for simulation and reconstruction
type ETree struct {
	Name         string
	BranchLength float64
	OnlyTopology bool
	Fasta        *fasta.Fasta //assigning fastas to nodes
	State        int //corresponds to a base same numbering encoded by dna.Base
	Stored       []float64 //a list of probabilities for each base at any given site of the genome
	Scrap        float64
	Left         *ETree
	Right        *ETree
	Up           *ETree //traversing the tree for reconstruction
}

//read tree from filename and add fastas and up pointers to the tree
func ReadTree(newickFilename string, fastasFilename string) *ETree {
	tr := ReadNewick(newickFilename)
	AssignFastas(tr, fastasFilename)
	return tr
}

//read in tree from filename
func ReadNewick(filename string) *ETree {
	// some small tweaks so that the code dies right away instead of continuing with error
	var line string

	file, err := os.Open(filename)
	common.ExitIfError(err)
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		if !strings.HasPrefix(line, "#") {
			return parseNewick(line[strings.Index(line, "("):(1 + strings.LastIndex(line, ";"))])
		}
	}
	log.Fatalf("Error: %s is either empty or has no non-comment lines\n", filename)
	return nil
}

// read in tree from string of newick
func NewickToTree(tree string) *ETree {
	makeTree := parseNewick(tree)
	return makeTree
}

func GetTree(node *ETree) []*ETree {
	var branch []*ETree
	branch = append(branch, node)
	if node.Right != nil {
		b := GetTree(node.Right)
		branch = append(branch, b...)
	}
	if node.Left != nil {
		a := GetTree(node.Left)
		branch = append(branch, a...)
	}
	return branch
}

func CopyTree(tree *ETree) *ETree {
	var treeCopy *ETree
	*treeCopy = *tree

	if tree.Left != nil {
		*treeCopy.Left = *CopyTree(tree.Left)
	}
	if tree.Right != nil {
		*treeCopy.Right = *CopyTree(tree.Right)
	}
	return treeCopy
}

func GetBranch(node *ETree) []*ETree {
	var branch []*ETree
	if node.Left != nil && node.Right != nil {
		branch = append(branch, node)
		a := GetBranch(node.Left)
		b := GetBranch(node.Right)
		branch = append(branch, a...)
		branch = append(branch, b...)
	}
	return branch
}

func GetLeaves(node *ETree) []*ETree {
	var leaf []*ETree
	if node.Left != nil && node.Right != nil {
		a := GetLeaves(node.Left)
		b := GetLeaves(node.Right)
		leaf = append(leaf, a...)
		leaf = append(leaf, b...)
	}
	if node.Left == nil && node.Right == nil {
		leaf = append(leaf, node)
	}
	return leaf
}

func splittingCommaIndex(input string) int {
	var openCount, closedCount = 0, 0
	for i, r := range input {
		if r == ',' && openCount == closedCount+1 {
			return i
		} else if r == '(' {
			openCount++
		} else if r == ')' {
			closedCount++
		}
	}
	return -1
}

func splitNameAndLength(input string) (string, float64, bool) {
	numColon := strings.Count(input, ":")
	if numColon == 0 {
		return input, 1, true
	} else if numColon == 1 {
		lastColon := strings.LastIndex(input, ":")
		name := input[:lastColon]
		branchLength := common.StringToFloat64(input[lastColon+1:])
		return name, branchLength, false
	}
	log.Fatalf("Error: %s should only have one or two colons\n", input)
	return "", 0, false
}

func parseNewickHelper(input string) *ETree {
	var answer ETree

	// all the character finding is up front and then all the logic follows
	// this is inefficient, but parsing trees should not be a limiting step
	// and I think this makes it more readable
	numOpen := strings.Count(input, "(")
	numClosed := strings.Count(input, ")")
	numComma := strings.Count(input, ",")
	numColon := strings.Count(input, ":")
	firstOpen := strings.Index(input, "(")
	lastClosed := strings.LastIndex(input, ")")
	lastColon := strings.LastIndex(input, ":")
	splittingComma := splittingCommaIndex(input)

	if numOpen != numClosed {
		log.Fatalf("Error: %s does not have an equal number of open and close parentheses\n", input)
	} else if numOpen != numComma {
		log.Fatalf("Error: %s does not have an a number of commas equal to the number of parenthesis pairs\n", input)
	} else if len(input) == 0 {
		log.Fatal("Error: can not build tree/node from an empty string")
	} else if numColon != 0 && (numColon != 2*numComma && lastColon < lastClosed) && (numColon != 2*numComma+1 && lastColon > lastClosed) {
		log.Fatalf("Error: %s should have a number of colons equal to zero or twice the number of colons (with another possible for the root branch)\n", input)
	}

	answer = ETree{Name: "", BranchLength: 1, OnlyTopology: true, Fasta: nil, State: 4, Stored: []float64{0, 0, 0, 0}, Scrap: 0.0, Left: nil, Right: nil, Up: nil}
	if numOpen == 0 { /* leaf node */
		answer.Name, answer.BranchLength, answer.OnlyTopology = splitNameAndLength(input)
		return &answer
	} else { /* internal node */
		answer.Name, answer.BranchLength, answer.OnlyTopology = splitNameAndLength(input[lastClosed+1:])
		answer.Left = parseNewickHelper(input[firstOpen+1 : splittingComma])
		answer.Right = parseNewickHelper(input[splittingComma+1 : lastClosed])
		return &answer
	}
}

func parseNewick(input string) *ETree {
	if !strings.HasPrefix(input, "(") || !strings.HasSuffix(input, ";") {
		log.Fatalf("Error: tree %s should start with '(' and end with ';'", input)
	}
	answer := parseNewickHelper(input[:len(input)-1])
	setUp(answer, nil)
	return answer
}

//tell tree what "up" is
func setUp(root *ETree, prevNode *ETree) {
	root.Up = prevNode
	if root.Left != nil && root.Right != nil {
		setUp(root.Left, root)
		setUp(root.Right, root)
	}
}

//set up tree with fastas
func AssignFastas(root *ETree, fastaFilename string) {
	fastas := fasta.Read(fastaFilename)
	//SetUp(root, nil) placed inside the parseNewick function so it always happens
	leaves := GetLeaves(root)
	for i := 0; i < len(leaves); i++ {
		for j := 0; j < len(fastas); j++ {
			if leaves[i].Name == fastas[j].Name {
				leaves[i].Fasta = fastas[j]
				leaves[i].State = int(leaves[i].Fasta.Seq[0])
			}
		}
		if leaves[i].Fasta == nil {
			log.Fatalf("Error: could not find %s in the fasta file\n", leaves[i].Name)
		}
	}
	branches := GetBranch(root)

	for i := 0; i < len(branches); i++ {
		f := fasta.Fasta{branches[i].Name, []dna.Base{}}
		branches[i].Fasta = &f
	}
}
