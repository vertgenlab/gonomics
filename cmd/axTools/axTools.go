// Command Group: "General Tools"

// Utilities for axt alignments
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"log"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func usage() {
	fmt.Print(
		"axTools - utilities for axt alignments\n" +
			"Usage:\n" +
			" axTools [options] input.axt output.axt/txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	var targetGaps *bool = flag.Bool("gap", false, "Find axt alignments when target contains Ns, but query does not")
	var querySwap *bool = flag.Bool("swap", false, "Swap target and query records. Must provide a `target.sizes and query.sizes` file containing query sequence lengths")

	var tLen *string = flag.String("tLen", "", "target `chrom.sizes` file containing target sequence lengths")
	var qLen *string = flag.String("qLen", "", "query `chrom.sizes` file containing query sequence lengths")

	var consensus *string = flag.String("fasta", "", "Output `.fa` consensus sequence based on the axt alignment")
	var minScore *int = flag.Int("minScore", 0, "filter axt alignments by minimum score")

	var stats *bool = flag.Bool("stats", false, "Calculate summary statistics (Alignment size, percent identity) for an AXT file. Can be used with the -bedfile option. Provide a .txt file in the output file argument that the stats will be sent to.")
	var bedfile *string = flag.String("bedfile", "", "For use with -stats. Provide a bed file, and only AXT records overlapping a bed record (relationship: any) will be used to calculate stats.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	input, output := flag.Arg(0), flag.Arg(1)

	if *bedfile != "" && !*stats {
		log.Fatalf("ERROR: -bedfile must be used with -stats.\n")
	}

	if *targetGaps {
		filterAxt(input, output)
	} else if fasta.IsFasta(*consensus) {
		axtToFa(input, output, *consensus)
	} else if *querySwap {
		QuerySwapAll(input, output, *tLen, *qLen)
	} else if *minScore != 0 {
		filterAxtScore(input, output, *minScore)
	} else if *stats {
		axtStats(input, output, *bedfile)
	} else {
		flag.Usage()
		if len(flag.Args()) != expectedNumArgs {
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
	}
}

// axtStats() takes in an AXT file, and creates an output text file with the length (of the reference interval) and percent ID
// (to the query). Optionally, a bed file can be provided, and only AXT records overlapping a bed record will be analyzed.
func axtStats(inputAxtFilename string, bedfile string, outputStatsFilename string) {
	var length int
	var pID float64
	var bedTree map[string]*interval.IntervalNode

	//read AXT to channel
	data, _ := axt.GoReadToChan(inputAxtFilename)
	out := fileio.EasyCreate(outputStatsFilename)
	fileio.WriteToFileHandle(out, "length\tpercentIdentity")

	//if a bed file is provided, read the bedfile and create a search tree for interval query
	if bedfile != "" {
		bd := bed.Read(bedfile)
		bedTree = interval.BedSliceToIntervalMap(bd)
	}

	//loop through AXT file
	for i := range data {
		if bedfile != "" {
			//if a bed file is provided, check to see if the AXT record overlaps a bed region
			if !interval.QueryBool(bedTree, i, "any", []interval.Interval{}) {
				//returns false if no overlap, so we want to continue
				continue
			}
		}
		//made it through checks, calculate stats for the AXT and write to file
		length, pID = calcAxtStat(i)
		fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%.2f", length, pID))
	}

	exception.PanicOnErr(out.Close())
}

// calcAxtStat() takes in an AXT record to calculate size and percent ID
func calcAxtStat(a axt.Axt) (int, float64) {
	return interval.IntervalSize(a), lift.AxtPercentIdentityInInterval(a, a)
}

func filterAxt(input string, output string) {
	ioWriter := fileio.EasyCreate(output)
	data, _ := axt.GoReadToChan(input)

	var index int = 0
	for each := range data {
		if axtTargetGap(each) {
			axt.WriteToFileHandle(ioWriter, each, index)
			index++
		}
	}
}

func axtTargetGap(record axt.Axt) bool {
	if dna.CountBase(record.RSeq, dna.N) != 0 && dna.CountBase(record.QSeq, dna.N) == 0 {
		return true
	} else {
		return false
	}
}

func axtQueryGap(record axt.Axt) bool {
	if dna.CountBase(record.QSeq, dna.N) != 0 && dna.CountBase(record.RSeq, dna.N) == 0 {
		return true
	} else {
		return false
	}
}

func axtToFa(input string, output string, target string) {
	ioWriter := fileio.EasyCreate(output)
	faMap := fasta.ToMap(fasta.Read(target))
	data, _ := axt.GoReadToChan(input)

	for each := range data {
		fasta.WriteFasta(ioWriter, axtSeq(each, faMap[each.RName]), 50)
	}
}

func filterAxtScore(input string, output string, minScore int) {
	ioWriter := fileio.EasyCreate(output)
	data, _ := axt.GoReadToChan(input)

	var index int
	for each := range data {
		if each.Score >= minScore {
			axt.WriteToFileHandle(ioWriter, each, index)
			index++
		}
	}
}

// if target sequence contains Ns, uses query non N bases to fill Ns.
func axtSeq(axtRecord axt.Axt, faSeq []dna.Base) fasta.Fasta {
	consensus := fasta.Fasta{
		Name: axtRecord.RName,
		Seq:  make([]dna.Base, 0, len(faSeq)),
	}
	consensus.Seq = append(consensus.Seq, faSeq[:axtRecord.RStart-1]...)
	for i := 0; i < len(axtRecord.RSeq); i++ {
		if axtRecord.RSeq[i] == dna.N && axtRecord.QSeq[i] != dna.N {
			consensus.Seq = append(consensus.Seq, axtRecord.QSeq[i])
		} else {
			consensus.Seq = append(consensus.Seq, axtRecord.RSeq[i])
		}
	}
	consensus.Seq = append(consensus.Seq, faSeq[axtRecord.REnd:]...)
	if len(consensus.Seq) != len(faSeq) {
		log.Fatalf("Error: Sequence length is not the same...\n")
	}
	return consensus
}

func QuerySwapAll(input string, output string, targetLen string, queryLen string) {
	targetInfo := chromInfo.ReadToMap(targetLen)
	queryInfo := chromInfo.ReadToMap(queryLen)

	axtWriter := fileio.EasyCreate(output)
	axtReader, _ := axt.GoReadToChan(input)

	var index int = 0
	for each := range axtReader {
		axt.Swap(&each, targetInfo[each.RName].Size, queryInfo[each.QName].Size)
		axt.WriteToFileHandle(axtWriter, each, index)
		index++
	}
	exception.PanicOnErr(axtWriter.Close())
}
