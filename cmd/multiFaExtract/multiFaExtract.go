// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func multiFaExtract(s Settings) {
	var ans []fasta.Fasta
	records := fasta.Read(s.InFile)
	if s.Bed == "" {
		ans = extractMultiHelper(records, s.Start, s.End, s.RemoveGaps, false)
		fasta.Write(s.OutFile, ans)
	} else {
		bedChan := bed.GoReadToChan(s.Bed)
		for b := range bedChan {
			ans = extractMultiHelper(records, b.ChromStart, b.ChromEnd, s.RemoveGaps, b.Strand == bed.Negative)
			if b.FieldsInitialized >= 4 {
				fasta.Write(fmt.Sprintf("%s.fa", b.Name), ans)
			} else {
				fasta.Write(fmt.Sprintf("%s.%d.%d.fa", b.Chrom, b.ChromStart, b.ChromEnd), ans)
			}
		}
	}
}

func extractMultiHelper(records []fasta.Fasta, start int, end int, removeGaps bool, revComp bool) []fasta.Fasta {
	var ans = make([]fasta.Fasta, len(records))
	if !(start < end) {
		log.Fatalf("Invalid arguments, start must be lower than end. start=%d end=%d\n", start, end)
	}
	for i := range records {
		ans[i] = fasta.Extract(records[i], fasta.RefPosToAlnPos(records[0], start), fasta.RefPosToAlnPos(records[0], end), records[i].Name)
		if revComp {
			ans[i] = fasta.Copy(ans[i])
			fasta.ReverseComplement(ans[i])
		}
	}
	if removeGaps {
		ans = fasta.RemoveGaps(ans)
	}
	return ans
}

func usage() {
	fmt.Print(
		"multiFaExtract - Pull sub-sequence from multiple Fasta alignment for each entry. Uses reference indices, treating the first sequence as the reference.\n" +
			"Usage:\n" +
			"multiFaExtract multi.fa out.fa start end\n" +
			"OR\n" +
			"multiFaExtract -bed regions.bed multi.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile     string
	OutFile    string
	Start      int
	End        int
	Bed        string
	RemoveGaps bool
}

func main() {
	var expectedNumArgs int = 4
	var removeGaps *bool = flag.Bool("removeGaps", false, "Removes gaps from the output sequences. Note that the output will no longer be in valid multiFa format.")
	var bed *string = flag.String("bed", "", "Extract all regions from start and end positions specified by the input bed file. Ignores the chrom field, so format the input bed file to contain only entries matching the multiFa chromosome. Output files will be named according to the coordinates of the bed entry.")
	var s Settings

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *bed != "" {
		expectedNumArgs = 1
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	if *bed == "" {
		s = Settings{
			InFile:     flag.Arg(0),
			OutFile:    flag.Arg(1),
			Start:      common.StringToInt(flag.Arg(2)),
			End:        common.StringToInt(flag.Arg(3)),
			Bed:        *bed,
			RemoveGaps: *removeGaps,
		}
	} else {
		s = Settings{
			InFile:     flag.Arg(0),
			OutFile:    "",
			Start:      -1,
			End:        -1,
			Bed:        *bed,
			RemoveGaps: *removeGaps,
		}
	}
	multiFaExtract(s)
}
