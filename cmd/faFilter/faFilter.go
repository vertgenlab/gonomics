package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faFilter(infile string, outfile string, name string, notName string, refPositions bool, start int, end int, minSize int) {
	records := fasta.Read(infile)
	var outlist []fasta.Fasta
	var pass bool = true

	if start > end && end != -1 {
		log.Fatalf("End must be larger than Start.")
	}

	if refPositions {
		start = fasta.RefPosToAlnPos(records[0], start)
		end = fasta.RefPosToAlnPos(records[0], end)
	}

	for i := 0; i < len(records); i++ {
		pass = true
		if name != "" && records[i].Name != name {
			pass = false
		}
		if notName != "" && records[i].Name == notName {
			pass = false
		}
		if len(records[i].Seq) < minSize {
			pass = false
		}
		if pass {
			if end == -1 {
				records[i].Seq = records[i].Seq[start:]
			} else {
				records[i].Seq = records[i].Seq[start:end]
			}
			outlist = append(outlist, records[i])
		}
	}
	fasta.Write(outfile, outlist)
}

func usage() {
	fmt.Print(
		"faFilter - Returns a filtered fasta based on argument parameters.\n" +
			"Usage:\n" +
			" faFindFast input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var refPositions *bool = flag.Bool("refPositions", false, "Uses reference positions for range specifications instead of alignment positions.")
	var start *int = flag.Int("start", 0, "Retains the sequence after this position.")
	var end *int = flag.Int("end", -1, "Retains the sequence before this position.")
	var name *string = flag.String("name", "", "Specifies the fasta record name.")
	var notName *string = flag.String("notName", "", "Returns all fasta records except for this input.")
	var minSize *int = flag.Int("minSize", 0, "Retains all fasta records with a sequence of at least that size")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFilter(inFile, outFile, *name, *notName, *refPositions, *start, *end, *minSize)
}
