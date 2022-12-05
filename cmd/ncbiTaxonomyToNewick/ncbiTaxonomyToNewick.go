// Command Group: "Sequence Evolution & Reconstruction"
// Command Usage: "Convert NCBI Taxonomy To Newick Format"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/tree/multifurtree"
	"log"
)

func prune(node *multifurtree.Node) *multifurtree.Node {
	var rank = map[string]int{
"superkingdom":1,
"kingdom":2,
"subkingdom":3,
"superphylum":4,
"phylum":5,
"subphylum":6,
"infraphylum":7,
"superclass":8,
"class":9,
"subclass":10,
"infraclass":11,
"cohort":12,
"subcohort":13,
"superorder":14,
"order":15,
"suborder":16,
"infraorder":17,
"parvorder":18,
"superfamily":19,
"family":20,
"subfamily":21,
"tribe":22,
"subtribe":23,
"genus":24,
"subgenus":25,
"section":26,
"subsection":27,
"series":28,
"subseries":29,
"species group":30,
"species subgroup":31,
"species":32,
"forma specialis":33,
"subspecies":34,
"varietas":35,
"subvariety":36,
"forma":37,
"serogroup":38,
"serotype":39,
"strain":40,
"isolate":41,
	}
	 minRank(
	return node
}

func minRank(node multifurtree.Node) int {
	var answer int = math.MaxInt32
	answer = numbers.Min(answer, node.BranchLength)
	for i := range node.Children {
		answer = numbers.Min(answer, minRank(node.Children[i]))
	}
	return answer
}



// ncbiTaxonomyToNewick takes in a nodes.dmp file and writes a newick file
func ncbiTaxonomyToNewick(nodesFilename string, newickFilename string) {
	var node multifurtree.Node
	var err error
	tree = multifurtree.ReadNewick(nodesFilename)
	err = multifurtree.WriteNewick(newickFilename, tree)
	if err != nil {
		log.Fatalf(err)
	}
}

func usage() {
	fmt.Print(
		"drawNewickTree takes in a newick format text file and outputs a png for tree visualization. If image width and height aren't specified they will default to 1500\n" +
			"Usage:\n" +
			"drawNewicktree [-option=int] <input.txt> <output.png>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2
	var imgWidth *int = flag.Int("imgWidth", 1500, "Specifies the width of the ouput.png.")
	var imgHeight *int = flag.Int("imgHeight", 1500, "Specifies the height of the ouput.png.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	imgOutFile := flag.Arg(1)

	drawNewickTree(infile, imgOutFile, *imgWidth, *imgHeight)
}
