package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"strings"
)

func crisprLib(inFile, outFile string, guidesPerOligo int) {
	var data bool = true
	var currGuide fasta.Fasta
	var currOligo []fasta.Fasta = []fasta.Fasta{{Name: ""}}

	upstream := dna.StringToBases("GGTCGAGCCGGAACTCGTCTCACACCG")
	downstream := dna.StringToBases("GTTTgGAGACGTCTGGGTGCGCATCC")
	spacer := dna.StringToBases("GTTTTGAGACGgactgcCGTCTCcCACCG")

	fa := fasta.GoReadToChan(inFile)
	out := fileio.EasyCreate(outFile)

	for data {
		var name []string
		currOligo[0] = fasta.Fasta{Name: "", Seq: upstream}
		for i := 0; i < guidesPerOligo; i++ {
			currGuide, data = <-fa
			if !data {
				break
			}
			name = append(name, currGuide.Name)
			if currGuide.Seq[0] == 2 {
				currGuide.Seq = currGuide.Seq[1:]
			}
			currOligo[0].Seq = append(currOligo[0].Seq, currGuide.Seq...)
			if i < guidesPerOligo-1 {
				currOligo[0].Seq = append(currOligo[0].Seq, spacer...)
			}
		}
		if !data {
			break
		}
		currOligo[0].Seq = append(currOligo[0].Seq, downstream...)
		currOligo[0].Name = strings.Join(name, "_")
		fasta.WriteToFileHandle(out, currOligo, numbers.MaxInt)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("\ncrisprLib -- format oligos for ordering for a library of CRISPR guides\n" +
		"Input a fasta file of guides (one guide per fasta entry). Spacers and BsmBI sites are added between guides that create overhangs suitable for golden gate cloning with most standard CRISPR vectors. " +
		"Primer sequences are added to the ends of oligos for amplification. It is recommended to use 2 or 3 guides per oligo.\n" +
		"Use these primers to amplify the oligo library: \n" +
		"Fwd: GGTCGAGCCGGAACTCGTC\n" +
		"Rev: GGATGCGCACCCAGACGTCTC\n\n" +
		"Usage:\n" +
		"\tcrisprLib [options] inGuides.fa oligosForOrder.fa\n\n")
}

func main() {
	var guidesPerOligo *int = flag.Int("guidesPerOligo", 2, "Number of CRISPR guides to be concatenated into an oligo for ordering. Default value is 2.")
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFa := flag.Arg(0)
	outOligos := flag.Arg(1)

	crisprLib(inFa, outOligos, *guidesPerOligo)
}
