// Command Group: "Sequence Evolution & Reconstruction"
// Command Usage: "Convert NCBI Taxonomy To Newick Format"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/tree/multifurtree"
	"log"
	"math"
	"strings"
)

var rank = map[string]float64{
	"superkingdom":     1,
	"kingdom":          2,
	"subkingdom":       3,
	"superphylum":      4,
	"phylum":           5,
	"subphylum":        6,
	"infraphylum":      7,
	"superclass":       8,
	"class":            9,
	"subclass":         10,
	"infraclass":       11,
	"cohort":           12,
	"subcohort":        13,
	"superorder":       14,
	"order":            15,
	"suborder":         16,
	"infraorder":       17,
	"parvorder":        18,
	"superfamily":      19,
	"family":           20,
	"subfamily":        21,
	"tribe":            22,
	"subtribe":         23,
	"genus":            24,
	"subgenus":         25,
	"section":          26,
	"subsection":       27,
	"series":           28,
	"subseries":        29,
	"species group":    30,
	"species subgroup": 31,
	"species":          32,
	"forma specialis":  33,
	"subspecies":       34,
	"varietas":         35,
	"subvariety":       36,
	"forma":            37,
	"serogroup":        38,
	"serotype":         39,
	"strain":           40,
	"isolate":          41,
}

func prune(node *multifurtree.Node, threshold string) *multifurtree.Node {
	//if node.BranchLength > rank[threshold] {
	//	return nil
	//}
	var newChildren []*multifurtree.Node = make([]*multifurtree.Node, 0)
	for i := range node.Children {
		if minRank(node.Children[i]) <= rank[threshold] {
			newChildren = append(newChildren, prune(node.Children[i], threshold))
		}
	}
	node.Children = newChildren
	return node
}

func minRank(node *multifurtree.Node) float64 {
	var answer float64 = math.MaxFloat32
	answer = numbers.Min(answer, node.BranchLength)
	for i := range node.Children {
		answer = numbers.Min(answer, minRank(node.Children[i]))
	}
	return answer
}

func readNodesDmp(filename string) *multifurtree.Node {
	var nodes map[string]*multifurtree.Node = make(map[string]*multifurtree.Node)
	var currRank float64
	var words []string
	var currId, parentId, line string
	var currNode, parentNode *multifurtree.Node
	var doneReading, found bool
	var err error
	file := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words = strings.Split(line, "\t")
		currId = words[0]
		parentId = words[2]
		currRank = rank[words[4]]
		currNode, found = nodes[currId]
		if !found {
			currNode = &multifurtree.Node{Name: currId, OnlyTopology: true, BranchLength: currRank, Children: nil}
			nodes[currId] = currNode
		} else {
			currNode.BranchLength = currRank
		}
		parentNode, found = nodes[parentId]
		if !found {
			parentNode = &multifurtree.Node{Name: parentId, OnlyTopology: true, Children: nil}
			nodes[parentId] = parentNode
		}
		if currId != parentId {
			parentNode.Children = append(parentNode.Children, currNode)
		}
	}
	err = file.Close()
	exception.PanicOnErr(err)
	fmt.Printf("%d\n", len(nodes["131567"].Children))
	return nodes["131567"] // return cellular organisms
}

// ncbiTaxonomyToNewick takes in a nodes.dmp file and writes a newick file
func ncbiTaxonomyToNewick(nodesFilename string, newickFilename string) {
	var tree *multifurtree.Node
	var err error
	tree = readNodesDmp(nodesFilename)
	tree = prune(tree, "superfamily")
	err = multifurtree.WriteNewick(newickFilename, tree)
	if err != nil {
		log.Fatalf("Error when writing newick file: %s\n", err)
	}
}

func usage() {
	fmt.Print(
		"ncbiTaxonomyToNewick takes a nodes.dmp file and returns a newick file\n" +
			"Usage:\n" +
			"drawNewicktree <nodes.dmp> <output.nh>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFilename := flag.Arg(0)
	outFilename := flag.Arg(1)

	ncbiTaxonomyToNewick(inFilename, outFilename)
}
