// Package expandedTree contains structs and functions that build and utilize information in a binary tree format that can be used for genome simulation and reconstruction
package expandedTree

import (
	"errors"
	"fmt"
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/tree"
)

// ETree is a struct that represents a node in a binary tree, and has additional fields for simulation and reconstruction.
type ETree struct {
	Name                  string
	BranchLength          float64 // This is the branch length between this node and Up.
	OnlyTopology          bool
	Fasta                 *fasta.Fasta // Contains the sequence associated with this node. Fasta.Name should equal Name.
	State                 int          // corresponds to a base
	Stored                []float64    // a list of probabilities for each base at any given site of the genome
	Scrap                 float64
	Left                  *ETree
	Right                 *ETree
	Up                    *ETree         // The immediate ancestral node.
	DescendentBasePresent bool        // True if any descendent nodes have a base, in a specific position
	BasePresent           bool        // True if this node has a base (A, C, G, T, or N). False if this node has dna.Gap.
	SubstitutionMatrix    [][]float64 // for custom substitution matrices. This is a 4x4 substitution matrices for nucleotides.
}

// ReadTree takes a filename of a tree in newick format and a filename of a fasta file.  The names in the tree will be assigned
// a sequence in the fasta file if they have identical names.  The function returns a pointer to the root node.
func ReadTree(newickFilename string, fastasFilename string) (*ETree, error) {
	tr, err := ReadNewick(newickFilename)
	if err != nil {
		return nil, err
	}
	AssignFastas(tr, fastasFilename)
	return tr, nil
}

// ReadNewick takes a filename of a tree in newick format and returns a pointer to the root node.
func ReadNewick(filename string) (*ETree, error) {
	var singleLineTree string
	singleLineTree = fileio.ReadFileToSingleLineString(filename)
	return parseNewick(singleLineTree[strings.Index(singleLineTree, "(") : 1+strings.LastIndex(singleLineTree, ";")])
}

// GetTree takes a root node and returns a slice of all nodes in the tree.
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

// CopyTree takes a root node and returns a pointer to the root node of a copy of the entire tree.
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

// GetBranch takes a root node and returns a slice of all internal nodes in the tree.
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

// GetLeaves takes a root node and returns a slice of all leaf nodes in the tree.
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

func splitNameAndLength(input string) (string, float64, bool, error) {
	numColon := strings.Count(input, ":")
	if numColon == 0 {
		return input, 1, true, nil
	} else if numColon == 1 {
		lastColon := strings.LastIndex(input, ":")
		name := input[:lastColon]
		branchLength, err := strconv.ParseFloat(input[lastColon+1:], 64)
		if err != nil {
			return name, 0, false, fmt.Errorf("Error: expecting %s to be a branch length\n", input[lastColon+1:])
		}
		return name, branchLength, false, nil
	}
	return "", 0, false, fmt.Errorf("Error: %s should only have one or two colons\n", input)
}

func parseNewickHelper(input string) (*ETree, error) {
	var answer ETree
	var err error

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
		return nil, fmt.Errorf("Error: %s does not have an equal number of open and close parentheses\n", input)
	} else if numOpen != numComma {
		return nil, fmt.Errorf("Error: %s does not have an a number of commas equal to the number of parenthesis pairs. Could be caused by a non-bifurcating tree\n", input)
	} else if len(input) == 0 {
		return nil, errors.New("Error: can not build tree/node from an empty string")
	} else if numColon != 0 && (numColon != 2*numComma && lastColon < lastClosed) && (numColon != 2*numComma+1 && lastColon > lastClosed) {
		return nil, fmt.Errorf("Error: %s should have a number of colons equal to zero or twice the number of colons (with another possible for the root branch)\n", input)
	}

	answer = ETree{Name: "", BranchLength: 1, OnlyTopology: true, Fasta: nil, State: 4, Stored: []float64{0, 0, 0, 0}, Scrap: 0.0, Left: nil, Right: nil, Up: nil}
	if numOpen == 0 { /* leaf node */
		answer.Name, answer.BranchLength, answer.OnlyTopology, err = splitNameAndLength(input)
		if err != nil {
			return nil, err
		}
		return &answer, nil
	} else { /* internal node */
		answer.Name, answer.BranchLength, answer.OnlyTopology, err = splitNameAndLength(input[lastClosed+1:])
		if err != nil {
			return nil, err
		}
		answer.Left, err = parseNewickHelper(input[firstOpen+1 : splittingComma])
		if err != nil {
			return nil, err
		}
		answer.Right, err = parseNewickHelper(input[splittingComma+1 : lastClosed])
		if err != nil {
			return nil, err
		}
		return &answer, nil
	}
}

func parseNewick(input string) (*ETree, error) {
	if !strings.HasPrefix(input, "(") || !strings.HasSuffix(input, ";") {
		return nil, fmt.Errorf("Error: tree %s should start with '(' and end with ';'", input)
	}
	return parseNewickHelper(input[:len(input)-1])
}

// SetUp takes a pointer to a node, root, and a "previous" node and sets the "Up" pointer
// of root to point to prevNode.
func SetUp(root *ETree, prevNode *ETree) {
	if prevNode != nil {
		root.Up = prevNode
	} else {
		root.Up = nil
	}
	if root.Left != nil && root.Right != nil {
		SetUp(root.Left, root)
		SetUp(root.Right, root)
	}
}

// AssignFastas takes a root node and the name of a fasta file.
// The function will assign fasta sequences to all the nodes
// in the tree if the names are identical.
func AssignFastas(root *ETree, fastaFilename string) {
	fastas := fasta.Read(fastaFilename)
	SetUp(root, nil)
	leaves := GetLeaves(root)
	for i := 0; i < len(leaves); i++ {
		for j := 0; j < len(fastas); j++ {
			if leaves[i].Name == fastas[j].Name {
				leaves[i].Fasta = &fastas[j]
				leaves[i].State = int(leaves[i].Fasta.Seq[0])
			}
		}
	}
	branches := GetBranch(root)

	for i := 0; i < len(branches); i++ {
		f := fasta.Fasta{branches[i].Name, []dna.Base{}}
		branches[i].Fasta = &f
	}
}

// FindNodeName takes a root node the name of the node we are searching for
// if the name is found at the root or below, a pointer to that node
// is the return value, and nil otherwise.
func FindNodeName(node *ETree, findMe string) *ETree {
	var searchResult *ETree
	if node == nil {
		return nil
	}
	if node.Name == findMe {
		return node
	}
	if node.Left != nil {
		searchResult = FindNodeName(node.Left, findMe)
		if searchResult != nil {
			return searchResult
		}
	}
	if node.Right != nil {
		searchResult = FindNodeName(node.Right, findMe)
		if searchResult != nil {
			return searchResult
		}
	}
	return nil
}

// ToNewickString converts a *ETree to a Newick-format string.
func ToNewickString(node *ETree) string {
	treeToWrite := toTree(node)
	return tree.ToString(treeToWrite)
}

// ToNewickFile writes a newick tree to a filename from a specified root.
func ToNewickFile(filename string, root *ETree) {
	out := fileio.EasyCreate(filename)
	_, err := fmt.Fprintf(out, "%v\n", ToNewickString(root))
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

// toTree converts an input *ETree to a *tree.Tree
func toTree(node *ETree) *tree.Tree {
	var answer *tree.Tree = &tree.Tree{
		Name:         node.Name,
		OnlyTopology: node.OnlyTopology,
		BranchLength: node.BranchLength,
	}
	if node.Left != nil {
		answer.Left = toTree(node.Left)
	}
	if node.Right != nil {
		answer.Right = toTree(node.Right)
	}
	return answer
}

// ToMap creates a map[string]*Etree, mapping each node's name to the *ETree struct.
func ToMap(root *ETree) map[string]*ETree {
	var answer = make(map[string]*ETree)
	toMapHelper(root, answer)
	return answer
}

// toMapHelper is a helper function of ToMap, and assists in the creation of a map[string]*ETree.
func toMapHelper(node *ETree, answer map[string]*ETree) {
	answer[node.Name] = node
	if node.Left == nil && node.Right != nil {
		log.Fatalf("Error: poorly formed binary tree. Node: %v has a nil left node and a non-nil right node.\n", node.Name)
	}
	if node.Left != nil && node.Right == nil {
		log.Fatalf("Error: poorly formed binary tree. Node: %v has a non-nil left node and a nil right node.\n", node.Name)
	}
	if node.Left != nil && node.Right != nil {
		toMapHelper(node.Left, answer)
		toMapHelper(node.Right, answer)
	}
}
