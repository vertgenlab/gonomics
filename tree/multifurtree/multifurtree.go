package tree

import (
	"bytes"
	"fmt"
	"os"
)

type Node struct {
	Name         string
	OnlyTopology bool
	BranchLength float64
	Children         []*Node
}

func toStringHelper(buffer *bytes.Buffer, node *Node) {
	if len(node.Children) == 0 {
		if node.OnlyTopology {
			buffer.WriteString(node.Name)
		} else {
			buffer.WriteString(fmt.Sprintf("%s:%f", node.Name, node.BranchLength))
		}
	} else {
		buffer.WriteRune('(')
		toStringHelper(buffer, node.Children[0])
		for i:=1 ; i<len(node.Children); i++ {
			buffer.WriteRune(',')
			toStringHelper(buffer, node.Children[i])
		}
		if node.OnlyTopology {
			buffer.WriteString(fmt.Sprintf(")%s", node.Name))
		} else {
			buffer.WriteString(fmt.Sprintf(")%s:%f", node.Name, node.BranchLength))
		}
	}
}

func ToString(node *Node) string {
	if node == nil {
		return ""
	}
	var buffer bytes.Buffer
	toStringHelper(&buffer, node)
	buffer.WriteRune(';')
	return buffer.String()
}

func WriteNewick(filename string, node *Node) error {
	file, err := os.Create(filename)
	defer file.Close()
	if err != nil {
		return err
	}
	fmt.Fprintf(file, "%s\n", ToString(node))
	return nil
}
