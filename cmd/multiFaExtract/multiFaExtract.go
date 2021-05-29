// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func multiFaExtract(infile string, outfile string, start int, end int, removeGaps bool) {
	if !(start < end) {
		log.Fatalf("Invalid arguments, start must be lower than end")
	}
	records := fasta.Read(infile)
	var ans []fasta.Fasta
	for i := 0; i < len(records); i++ {
		ans = append(ans, fasta.Extract(records[i], fasta.RefPosToAlnPos(records[0], start), fasta.RefPosToAlnPos(records[0], end), records[i].Name))
	}
	if removeGaps {
		ans = fasta.RemoveGaps(ans)
	}

	fasta.Write(outfile, ans)
}

func usage() {
	fmt.Print(
		"multiFaExtract - Pull sub-sequence from multiple Fasta alignment for each entry. Uses reference indices, treating the first sequence as the reference.\n" +
			"Usage:\n" +
			"multiFaExtract multi.fa out.fa start end\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var removeGaps *bool = flag.Bool("removeGaps", false, "Removes gaps from the output sequences. Note that the output will no longer be in valid multiFa format.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	start := common.StringToInt(flag.Arg(2))
	end := common.StringToInt(flag.Arg(3))

	multiFaExtract(infile, outfile, start, end, *removeGaps)
}
