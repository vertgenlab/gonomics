package tree

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Tree struct {
	Name         string
	OnlyTopology bool
	BranchLength float64
	Left         *Tree
	Right        *Tree
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

func parseNewick(input string) (*Tree, error) {
	if !strings.HasPrefix(input, "(") || !strings.HasSuffix(input, ";") {
		return nil, fmt.Errorf("Error: tree %s should start with '(' and end with ';'", input)
	}
	return parseNewickHelper(input[:len(input)-1])
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

func parseNewickHelper(input string) (*Tree, error) {
	var answer Tree
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

	answer = Tree{Name: "", OnlyTopology: true, BranchLength: 1, Left: nil, Right: nil}
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

func ReadNewick(filename string) (*Tree, error) {
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
			return parseNewick(line[strings.Index(line, "("):strings.LastIndex(line, ";")])
		}
	}
	return nil, errors.New("Error: tree file is either empty or has no non-comment lines")
}

func toStringHelper(buffer *bytes.Buffer, node *Tree) {
	if node.Left == nil && node.Right == nil {
		if node.OnlyTopology {
			buffer.WriteString(node.Name)
		} else {
			buffer.WriteString(fmt.Sprintf("%s:%f", node.Name, node.BranchLength))
		}
	} else {
		buffer.WriteRune('(')
		if node.Left != nil {
			toStringHelper(buffer, node.Left)
		}
		buffer.WriteRune(',')
		if node.Right != nil {
			toStringHelper(buffer, node.Right)
		}
		if node.OnlyTopology {
			buffer.WriteString(fmt.Sprintf(")%s", node.Name))
		} else {
			buffer.WriteString(fmt.Sprintf(")%s:%f", node.Name, node.BranchLength))
		}
	}
}

func ToString(node *Tree) string {
	if node == nil {
		return ""
	}
	var buffer bytes.Buffer
	toStringHelper(&buffer, node)
	buffer.WriteRune(';')
	return buffer.String()
}

func WriteNewick(filename string, node *Tree) error {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		return err
	}
	fmt.Fprintf(file, "%s\n", ToString(node))
	return nil
}
