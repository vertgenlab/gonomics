package main

import (
	"fmt"
	"log"
	"flag"
	"github.com/vertgenlab/gonomics/fasta"
)

func faFilter(infile string, outfile string, name string, refPositions bool, Start int64, End int64) {
	records := fasta.Read(infile)
	var outlist []*fasta.Fasta
	var pass bool = true

	if Start > End  && End != -1 {
		log.Fatalf("End must be larger than Start.")
	}

	if refPositions {
		Start = fasta.RefPosToAlnPos(records[0], Start)
		End = fasta.RefPosToAlnPos(records[0], End)
	}

	for i := 0; i < len(records); i++ {
		pass = true
		if name != "" {
			if records[i].Name != name {
				pass = false
			} else if End == -1 {
				records[i].Seq = records[i].Seq[Start:]
			} else {
				records[i].Seq = records[i].Seq[Start:End]
			}
		} else if End == -1 {
			records[i].Seq = records[i].Seq[Start:]
		} else {
			records[i].Seq = records[i].Seq[Start:End]
		}
		if pass {
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
	var Start *int64 = flag.Int64("start", 0, "Retains the sequence after this position.")
	var End *int64 = flag.Int64("end", -1, "Retains the sequence before this position.")
	var name *string = flag.String("name", "", "Specifies the fasta record name.")
	
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

	faFilter(inFile, outFile, *name, *refPositions, *Start, *End)
}