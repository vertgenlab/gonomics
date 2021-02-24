package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
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

	var concensus *string = flag.String("fasta", "", "Output `.fa` consensus sequence based on the axt alignment")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	input, output := flag.Arg(0), flag.Arg(1)

	if *targetGaps {
		filterAxt(input, output)
	} else if fasta.IsFasta(*concensus) {
		axtToFa(input, output, *concensus)
	} else if *querySwap {
		QuerySwapAll(input, output, *tLen, *qLen)
	} else {
		flag.Usage()
		if len(flag.Args()) != expectedNumArgs {
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
	}
}

func filterAxt(input string, output string) {
	ioWriter := fileio.EasyCreate(output)
	data := axt.GoReadToChan(input)

	var index int = 0
	for each := range data {
		if axtTargetGap(each) {
			axt.WriteToFileHandle(ioWriter, each, index)
			index++
		}
	}
}

func axtTargetGap(record *axt.Axt) bool {
	if dna.CountBase(record.RSeq, dna.N) != 0 && dna.CountBase(record.QSeq, dna.N) == 0 {
		return true
	} else {
		return false
	}
}

func axtQueryGap(record *axt.Axt) bool {
	if dna.CountBase(record.QSeq, dna.N) != 0 && dna.CountBase(record.RSeq, dna.N) == 0 {
		return true
	} else {
		return false
	}
}

func axtToFa(input string, output string, target string) {
	ioWriter := fileio.EasyCreate(output)
	faMap := fasta.ToMap(fasta.Read(target))
	data := axt.GoReadToChan(input)

	for each := range data {
		fasta.WriteFasta(ioWriter, axtSeq(each, faMap[each.RName]), 50)
	}
}

//if target sequence contains Ns, uses query non N bases to fill Ns
func axtSeq(axtRecord *axt.Axt, faSeq []dna.Base) *fasta.Fasta {
	concensus := &fasta.Fasta{
		Name: fmt.Sprintf("%s", axtRecord.RName),
		Seq:  make([]dna.Base, 0, len(faSeq)),
	}
	concensus.Seq = append(concensus.Seq, faSeq[:axtRecord.RStart-1]...)
	for i := 0; i < len(axtRecord.RSeq); i++ {
		if axtRecord.RSeq[i] == dna.N && axtRecord.QSeq[i] != dna.N {
			concensus.Seq = append(concensus.Seq, axtRecord.QSeq[i])
		} else {
			concensus.Seq = append(concensus.Seq, axtRecord.RSeq[i])
		}
	}
	concensus.Seq = append(concensus.Seq, faSeq[axtRecord.REnd:]...)
	if len(concensus.Seq) != len(faSeq) {
		log.Fatalf("Error: Sequence length is not the same...\n")
	}
	return concensus
}

func QuerySwapAll(input string, output string, targetLen string, queryLen string) {
	targetInfo := chromInfo.ReadToMap(targetLen)
	queryInfo := chromInfo.ReadToMap(queryLen)

	axtWriter := fileio.EasyCreate(output)
	axtReader := axt.GoReadToChan(input)

	var index int = 0
	var curr *axt.Axt
	for each := range axtReader {
		curr = axt.SwapBoth(each, int(targetInfo[each.RName].Size), int(queryInfo[each.QName].Size))
		axt.WriteToFileHandle(axtWriter, curr, index)
		index++
	}
	axtWriter.Close()
}
