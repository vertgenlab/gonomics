package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type restrictionEnzyme struct {
	name       string
	cutString  string
	cutBases   []dna.Base
	cutBasesRC []dna.Base
	cutPos     int
	pal        bool //palindrome
}

func getCutPos(cutSite string) (cutPos int, cutString []string) {
	var found bool = false
	for i := range cutSite {
		if string(cutSite[i]) == "^" {
			found = true
			cutPos = i
		} else {
			cutString = append(cutString, string(cutSite[i]))
		}
	}
	if found == false {
		log.Fatalf("The input restriction enzyme cut site must have the '^' character to denote the cut location. Your seq: %s", cutSite)
	}
	return cutPos, cutString
}

func digestGenome(genome string, cutSite string, outFile string) {
	var currChrom string
	var cutString []string
	var base, prevCut, numCut int
	var bedRegion bed.Bed
	var re restrictionEnzyme

	out := fileio.EasyCreate(outFile)

	//default enzymes
	switch cutSite {
	case "MboI":
		re.name = cutSite
		re.cutPos = 0
		re.cutBases = []dna.Base{2, 0, 3, 1}
		re.cutBasesRC = dna.ReturnReverseComplement(re.cutBases)
		re.pal = true
	case "DnpII":
		re.name = cutSite
		re.cutPos = 0
		re.cutBases = []dna.Base{2, 0, 3, 1}
		re.cutBasesRC = dna.ReturnReverseComplement(re.cutBases)
		re.pal = true
	case "BglII":
		re.name = cutSite
		re.cutPos = 1
		re.cutBases = []dna.Base{0, 2, 0, 3, 1, 3}
		re.cutBasesRC = dna.ReturnReverseComplement(re.cutBases)
		re.pal = true
	case "HindIII":
		re.name = cutSite
		re.cutPos = 1
		re.cutBases = []dna.Base{0, 0, 2, 1, 3, 3}
		re.cutBasesRC = dna.ReturnReverseComplement(re.cutBases)
		re.pal = true
	default:
		re.name = cutSite
		re.cutPos, cutString = getCutPos(cutSite)
		re.cutBases = dna.StringToBases(strings.Join(cutString, ""))
		re.cutBasesRC = dna.ReturnReverseComplement(re.cutBases)
		if dna.CompareSeqsIgnoreCase(re.cutBases, re.cutBasesRC) == 0 {
			re.pal = true
		} else {
			re.pal = false
		}
	}

	if re.pal == false {
		fmt.Println("Your input restriction enzyme cutting motif is not palindromic. Non-palindromtic motifs are currently not supported. Only the forward motif sequence will be searched.")
	}

	fa := fasta.GoReadToChan(genome)
	for i := range fa {
		prevCut = 0
		numCut = 0
		currChrom = i.Name
		for base = range i.Seq {
			testBases := i.Seq[base : base+len(re.cutBases)]
			if dna.CompareSeqsIgnoreCase(testBases, re.cutBases) == 0 {
				//found a cut site, write to bed
				bedRegion = bed.Bed{Chrom: currChrom, ChromStart: prevCut, ChromEnd: base + re.cutPos, Name: fmt.Sprintf("%s_%s_%d", cutSite, currChrom, numCut), FieldsInitialized: 4}
				bed.WriteToFileHandle(out, bedRegion)
				prevCut = base + re.cutPos
				numCut++
			}
		}
		//write the last fragment
		bedRegion = bed.Bed{Chrom: currChrom, ChromStart: prevCut, ChromEnd: len(i.Seq), Name: fmt.Sprintf("%s_%s_%d", cutSite, currChrom, numCut), FieldsInitialized: 4} //is len(i.seq) going to take too long/energy expensive. would love to do this without adding a chromSizes file
		bed.WriteToFileHandle(out, bedRegion)
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {

}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	var expectedNumArgs int = 3

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	digestGenome(flag.Arg(0), flag.Arg(1), flag.Arg(2))
}
