package tree_newick

import (
	"fmt"

	"bufio"
	"errors"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"strconv"
	"strings"
)

//Tree structure for simulation and reconstruction
type NTree struct {
	Name         string
	BranchLength float64
	OnlyTopology bool
	Fasta        *fasta.Fasta
	State        int
	Stored       []float64
	Scrap        float64
	Left         *NTree
	Right        *NTree
	Up           *NTree
}

//read tree from filename and add fastas and up pointers to the tree
func Read_tree(filename_newick string, filename_fastas string) *NTree {
	tr, err := ReadNewick(filename_newick)
	if err == nil {
	}
	Set_fastas_up(tr, filename_fastas)
	return tr
}

//read in tree from filename
func ReadNewick(filename string) (*NTree, error) {
	var line string

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		if !strings.HasPrefix(line, "#") {
			return parseNewick(line[strings.Index(line, "("):(1 + strings.LastIndex(line, ";"))])
		}
	}
	return nil, errors.New("Error: tree file is either empty or has no non-comment lines")
}

// read in tree from string of newick
func ReadNewick_string(tree string) *NTree {
	tree_thing, err := parseNewick(tree)
	if err != nil {
	}
	return tree_thing
}

func Get_tree(node *NTree) []*NTree {
	var branch []*NTree
	branch = append(branch, node)
	if node.Right != nil {
		b := Get_tree(node.Right)
		branch = append(branch, b...)
	}
	if node.Left != nil {
		a := Get_tree(node.Left)
		branch = append(branch, a...)
	}
	return branch
}

func Get_branch(node *NTree) []*NTree {
	var branch []*NTree
	if node.Left != nil && node.Right != nil {
		branch = append(branch, node)
		a := Get_branch(node.Left)
		b := Get_branch(node.Right)
		branch = append(branch, a...)
		branch = append(branch, b...)
	}
	return branch
}

func Get_leaf(node *NTree) []*NTree {
	var leaf []*NTree
	if node.Left != nil && node.Right != nil {
		a := Get_leaf(node.Left)
		b := Get_leaf(node.Right)
		leaf = append(leaf, a...)
		leaf = append(leaf, b...)
	}
	if node.Left == nil && node.Right == nil {
		leaf = append(leaf, node)
	}
	return leaf
}

func splittingCommaIndex(input string) int {
	var openCount, closedCount int = 0, 0
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

func parseNewickHelper(input string) (*NTree, error) {
	var answer NTree
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
		return nil, fmt.Errorf("Error: %s does not have an a number of commas equal to the number of parenthesis pairs\n", input)
	} else if len(input) == 0 {
		return nil, errors.New("Error: can not build tree/node from an empty string")
	} else if numColon != 0 && (numColon != 2*numComma && lastColon < lastClosed) && (numColon != 2*numComma+1 && lastColon > lastClosed) {
		return nil, fmt.Errorf("Error: %s should have a number of colons equal to zero or twice the number of colons (with another possible for the root branch)\n", input)
	}

	answer = NTree{Name: "", BranchLength: 1, OnlyTopology: true, Fasta: nil, State: 4, Stored: []float64{0, 0, 0, 0}, Scrap: 0.0, Left: nil, Right: nil, Up: nil}
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

func parseNewick(input string) (*NTree, error) {
	if !strings.HasPrefix(input, "(") || !strings.HasSuffix(input, ";") {
		return nil, fmt.Errorf("Error: tree %s should start with '(' and end with ';'", input)
	}
	return parseNewickHelper(input[:len(input)-1])
}

//set up tree with ups
func Set_up(root *NTree, prev_node *NTree) {
	if prev_node != nil {
		root.Up = prev_node
	} else {
		root.Up = nil
	}
	if root.Left != nil && root.Right != nil {
		Set_up(root.Left, root)
		Set_up(root.Right, root)
	}
}

//set up tree with fastas
func Set_fastas_up(root *NTree, fasta_file string) {
	fastas, err := fasta.Read(fasta_file)
	Set_up(root, nil)
	if err != nil {
	}
	leaves := Get_leaf(root)
	for i := 0; i < len(leaves); i++ {
		for j := 0; j < len(fastas); j++ {
			if leaves[i].Name == fastas[j].Name {
				leaves[i].Fasta = &fastas[j]
				leaves[i].State = int(leaves[i].Fasta.Seq[0])
			}

		}
	}
	branches := Get_branch(root)

	for i := 0; i < len(branches); i++ {
		f := fasta.Fasta{branches[i].Name, []dna.Base{}}
		branches[i].Fasta = &f
	}
}
