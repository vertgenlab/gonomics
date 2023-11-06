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
	var revMatch bool = false
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
		re.cutBasesRC = dna.ReverseComplementAndCopy(re.cutBases)
		re.pal = true
	case "DnpII":
		re.name = cutSite
		re.cutPos = 0
		re.cutBases = []dna.Base{2, 0, 3, 1}
		re.cutBasesRC = dna.ReverseComplementAndCopy(re.cutBases)
		re.pal = true
	case "BglII":
		re.name = cutSite
		re.cutPos = 1
		re.cutBases = []dna.Base{0, 2, 0, 3, 1, 3}
		re.cutBasesRC = dna.ReverseComplementAndCopy(re.cutBases)
		re.pal = true
	case "HindIII":
		re.name = cutSite
		re.cutPos = 1
		re.cutBases = []dna.Base{0, 0, 2, 1, 3, 3}
		re.cutBasesRC = dna.ReverseComplementAndCopy(re.cutBases)
		re.pal = true
	default:
		re.name = cutSite
		re.cutPos, cutString = getCutPos(cutSite)
		re.cutBases = dna.StringToBases(strings.Join(cutString, ""))
		re.cutBasesRC = dna.ReverseComplementAndCopy(re.cutBases)
		if dna.CompareSeqsIgnoreCase(re.cutBases, re.cutBasesRC) == 0 {
			re.pal = true
		} else {
			re.pal = false
		}
	}

	fa := fasta.GoReadToChan(genome)
	for i := range fa {
		prevCut = 0 //set the first fragment start base to 0
		numCut = 0
		currChrom = i.Name
		for base = range i.Seq { //loop over every base in the fasta file
			testBases := i.Seq[base : base+len(re.cutBases)] //get the next base in the fasta file plus the number of bases in the recognition motif. I don't get an index out of range error at the end of the sequence--not sure why
			if re.pal == false {
				if dna.CompareSeqsIgnoreCase(testBases, re.cutBasesRC) == 0 {
					revMatch = true
				}
			}
			if dna.CompareSeqsIgnoreCase(testBases, re.cutBases) == 0 || base+1 == len(i.Seq) || revMatch == true {
				//found a cut site or last fragment, write to bed
				switch {
				case base+1 == len(i.Seq):
					bedRegion = bed.Bed{Chrom: currChrom, ChromStart: prevCut, ChromEnd: len(i.Seq), Name: fmt.Sprintf("%s_%s_%d", cutSite, currChrom, numCut), Strand: '+', FieldsInitialized: 6}
				case revMatch == true:
					bedRegion = bed.Bed{Chrom: currChrom, ChromStart: prevCut, ChromEnd: base + (len(re.cutBases) - re.cutPos), Name: fmt.Sprintf("%s_%s_%d", cutSite, currChrom, numCut), Strand: '-', FieldsInitialized: 6}
					prevCut = base + (len(re.cutBases) - re.cutPos) //set the start of the next bed region to the end of the current bed region
					revMatch = false
				default:
					bedRegion = bed.Bed{Chrom: currChrom, ChromStart: prevCut, ChromEnd: base + re.cutPos, Name: fmt.Sprintf("%s_%s_%d", cutSite, currChrom, numCut), Strand: '+', FieldsInitialized: 6}
					prevCut = base + re.cutPos //set the start of the next bed region to the end of the current bed region
				}
				bed.WriteToFileHandle(out, bedRegion)
				numCut++
			}
		}
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("digestGenome -- Create a bed file of restriction fragments from an input FASTA file and a restriction enzyme recognition sequence\n" +
		"Input restriction enzyme recognition sequence must have the '^' to specify where the cut occurs. Example: CA^GT\n" +
		"The strand corresponds to whether the cut site was found on the forward or reverse strand. In the case of palindromic motifs, the positive strand and cut position will always be used." +
		"Several enzymes commonly used in Hi-C library preps are provided as defaults, and their name can be provided instead of the sequence motif:\n" +
		"MboI : ^GATC\n" +
		"DpnII : ^GATC\n" +
		"BglII : A^GATCT\n" +
		"HindIII : A^AGCTT\n" +
		"Example command with a default enzyme: digestGenome in.fa MboI out.bed\n" +
		"Usage:\n" +
		"digestGenome in.fa motif/default out.bed\n")
	flag.PrintDefaults()
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
