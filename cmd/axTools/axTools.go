// Command Group: "General Tools"

package main

import (
	"flag"
	"fmt"
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
			" axTools [options] input.axt output.axt\n" +
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
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	input, output := flag.Arg(0), flag.Arg(1)

	if *targetGaps {
		filterAxt(input, output)
	} else if fasta.IsFasta(*consensus) {
		axtToFa(input, output, *consensus)
	} else if *querySwap {
		QuerySwapAll(input, output, *tLen, *qLen)
	} else if *minScore != 0 {
		filterAxtScore(input, output, *minScore)
	} else {
		flag.Usage()
		if len(flag.Args()) != expectedNumArgs {
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
	}
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

// if target sequence contains Ns, uses query non N bases to fill Ns
func axtSeq(axtRecord axt.Axt, faSeq []dna.Base) fasta.Fasta {
	consensus := fasta.Fasta{
		Name: fmt.Sprintf("%s", axtRecord.RName),
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
	axtWriter.Close()
}
