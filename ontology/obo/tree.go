package obo

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

// BuildTree takes a slice of Obo structs and builds a tree, setting
// the Parents and Children fields to point to related nodes.
func BuildTree(terms []Obo) map[string]*Obo {
	termMap := makeTermMap(terms)

	var parentId string
	var foundInMap bool
	var parentTerm *Obo
	var isAItem IsADescription

	// Now we iterate through terms and set up the tree
	for i := range terms {
		term := &terms[i] //get a memory copy
		for _, isAItem = range term.IsA {
			parentId = isAItem.ParentId
			//fmt.Printf("parentId: %v\n", parentId)
			if parentTerm, foundInMap = termMap[parentId]; foundInMap {
				term.Parents = append(term.Parents, parentTerm)
				parentTerm.Children = append(parentTerm.Children, term)
			} else {
				//log.Fatalf("Error: The term with ID \"%s\" has a parent with ID \"%s\" that "+
				//	"is not found in the Obo file.\n", term.Id, parentId)
			}
		}
	}
	return termMap
}

// termToDot is a helper function of ToDot that converts and individual Obo struct to a
// DOT format digraph descriptor, and recursively does so for all children of the input Obo struct.
func termToDot(term Obo, out *fileio.EasyWriter, visitedNode map[string]bool) {
	var err error
	if visitedNode[term.Id] {
		return
	}
	visitedNode[term.Id] = true
	_, err = fmt.Fprintf(out, "\"%s\" [label = \"%s\"];\n", term.Id, term.Name)
	exception.PanicOnErr(err)

	// now we recursively visit children
	for _, child := range term.Children {
		_, err = fmt.Fprintf(out, "\"%s\" -> \"%s\";\n", term.Id, child.Id)
		termToDot(*child, out, visitedNode)
	}
}

// ToDot formats an input Obo slice as a tree and writes the resulting tree to a file in DOT format.
func ToDot(outFile string, terms []Obo) {
	var err error
	out := fileio.EasyCreate(outFile)
	BuildTree(terms)

	//write header line for Dot
	_, err = fmt.Fprintf(out, "digraph G{\n")
	exception.PanicOnErr(err)

	var visitedNode = make(map[string]bool)

	for _, term := range terms {
		termToDot(term, out, visitedNode)
	}

	//write tail line for Dot and close file
	_, err = fmt.Fprintf(out, "}\n")
	err = out.Close()
	exception.PanicOnErr(err)
}

// makeTermMap creates a map to find nodes based on their ID and AltIDs
func makeTermMap(terms []Obo) map[string]*Obo {
	var altId string
	termMap := make(map[string]*Obo)
	for i := range terms {
		term := &terms[i] // get a memory copy
		if termMap[term.Id] != nil {
			log.Fatalf("Term \"%v\" is duplicated in input slice of Obo.\n", term.Id)
		}
		termMap[term.Id] = term
		for _, altId = range term.AltIds {
			if termMap[altId] != nil {
				log.Fatalf("Term \"%v\" is duplicated in input slice of Obo.\n", term.Id)
			}
			termMap[altId] = term
		}
	}

	return termMap
}

// SubtreeToDot makes a DOT format tree for the subtree rooted at the input nodeId.
// BuildTree must have been run on this tree already, as termMap is an argument.
func SubtreeToDot(outFile string, nodeId string, termMap map[string]*Obo) {
	var err error
	out := fileio.EasyCreate(outFile)
	//write header line for Dot
	_, err = fmt.Fprintf(out, "digraph G{\n")
	exception.PanicOnErr(err)

	var visitedNode = make(map[string]bool)
	if _, ok := termMap[nodeId]; !ok {
		log.Fatalf("Input term not found in Obo file. Term: %v\n", nodeId)
	}
	termToDot(*termMap[nodeId], out, visitedNode)

	//write tail line for Dot and close file
	_, err = fmt.Fprintf(out, "}\n")
	err = out.Close()
	exception.PanicOnErr(err)
}

// NumberOfDescendents edits the field SubTreeSize of input Obo structs
// to the number of descendent nodes in a subtree rooted on a particular node.
func NumberOfDescendents(termMap map[string]*Obo) {
	var visitedNode = make(map[string]bool)
	for _, i := range termMap {
		numberOfDescendentsRecursive(i, visitedNode)
	}
}

// numberOfDescendentsRecursive is a helper function of NumberOfDescendents, which
// calculates the size of subtrees rooted on a particular Obo struct in a larger
// Obo tree.
func numberOfDescendentsRecursive(term *Obo, visitedNode map[string]bool) {
	if _, foundInMap := visitedNode[term.Id]; foundInMap {
		return
	}
	visitedNode[term.Id] = true
	//leaf case
	if len(term.Children) == 0 {
		term.SubTreeSize = 0
		return
	}
	for i := range term.Children {
		numberOfDescendentsRecursive(term.Children[i], visitedNode)
		//subtree size is the number of children plus the size of each child's subtree
		term.SubTreeSize += 1 + term.Children[i].SubTreeSize
	}
}

// SubTreeReport writes out the number of descendent nodes for each node in an Obo tree.
// This is used for debugging, and for picking subtrees to visualize based on size.
func SubTreeReport(outFile string, records []Obo) {
	var err error
	out := fileio.EasyCreate(outFile)
	for i := range records {
		_, err = fmt.Fprintf(out, "%v\tId: %v. Name: %v. Descendents: %v\n", records[i].SubTreeSize, records[i].Id, records[i].Name, records[i].SubTreeSize)
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
