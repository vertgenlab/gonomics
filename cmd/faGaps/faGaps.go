// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faNoGap(inFile string, outFile string) {
	records := fasta.Read(inFile)
	beds := bed.UngappedRegionsAllFromFa(records)
	bed.Write(outFile, beds)
}

func faSplitByNs(filename string, outFile string) {
	reader := fileio.EasyOpen(filename)
	writer := fileio.EasyCreate(outFile)
	defer reader.Close()
	defer writer.Close()
	for fa, done := fasta.NextFasta(reader); !done; fa, done = fasta.NextFasta(reader) {
		fasta.WriteToFileHandle(writer, chrSplitByNs(fa), 50)
	}
}

func faNoGapThreshold(inFile string, outFile string, windowSize int, ungappedBaseThreshold int) {
	reader := fileio.EasyOpen(inFile)
	writer := fileio.EasyCreate(outFile)
	var currBatch []bed.Bed //UngappedRegionsThresholdFromFa returns all ungapped beds for each chromosome, so we write them in batches.

	for fa, done := fasta.NextFasta(reader); !done; fa, done = fasta.NextFasta(reader) {
		currBatch = bed.UngappedRegionsThresholdFromFa(fa, windowSize, ungappedBaseThreshold)
		for i := range currBatch {
			bed.WriteBed(writer, currBatch[i])
		}
	}

	reader.Close()
	writer.Close()
}

func usage() {
	fmt.Print(
		"faGap - a program to investigate regions containing gaps\n\n" +
			"Usage:\n" +
			"  faGap [options] in.fa out.file\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var bedThreshold = flag.Bool("noGapBedThreshold", false, "Find genomic coordinates containing regions outside gaps in bed format. Scans the genome in sliding windows. Windows with a number of unmapped bases above the threshold contribute to an ungapped bed element.")
	var bedOut = flag.Bool("noGapBed", false, "find genomic coordinates containing regions outside gaps `noGap.bed`")
	var faOut = flag.Bool("fasta", false, "split fasta into several records using gapped coordinates `.fa/.gz`")
	var windowSize = flag.Int("windowSize", 500, "Set the size of the sliding window for the noGapBedThreshold function.")
	var ungappedBaseThreshold = flag.Int("ungappedBaseThreshold", 30, "Sets the number of ungapped bases that must be observed in a window for that window to contribute to an ungapped bed element.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	var inFile string = flag.Arg(0)
	var outFile string = flag.Arg(1)

	if *bedOut && *bedThreshold {
		log.Fatalf("Error: Cannot select both bedOut and bedThreshold. Use one or the other.")
	}

	if *bedThreshold {
		faNoGapThreshold(inFile, outFile, *windowSize, *ungappedBaseThreshold)
	} else if *bedOut {
		faNoGap(inFile, outFile)
	} else if *faOut {
		faSplitByNs(inFile, outFile)
	} else {
		log.Fatalf("Error: no subcommand specified.")
	}
}

//TODO: Belongs in the fasta package, but bed and fasta becomes an illegal import cycle
func chrSplitByNs(chr fasta.Fasta) []fasta.Fasta {
	unGapped := bed.UngappedRegionsFromFa(chr)
	var answer []fasta.Fasta = make([]fasta.Fasta, len(unGapped))
	for i := 0; i < len(unGapped); i++ {
		answer[i] = fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d", unGapped[i].Chrom, unGapped[i].ChromStart, unGapped[i].ChromEnd), Seq: chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]}
	}
	return answer
}
