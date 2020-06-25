package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
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

	var targetGaps *bool = flag.Bool("tGap", false, "Return '`axt` alignments when target contains Ns, but query does not")
	var queryGaps *bool = flag.Bool("qGap", false, "Return `.axt` alignments when query contains Ns, but target does not")
	//var concensus *bool = flag.Bool("concensus", false, "Output `.fasta` consensus sequence based on the axt alignment")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	input, output := flag.Arg(0), flag.Arg(1)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		//} else if *concensus {
		//	axtToFa(input, output)
	} else {
		searchAxt(input, output, targetGaps, queryGaps)
	}
}

func searchAxt(input string, output string, target *bool, query *bool) {
	file := fileio.EasyOpen(input)
	reader := make(chan *axt.Axt)
	writer := fileio.EasyCreate(output)
	go axt.ReadToChan(file, reader)
	var index int = 0
	for each := range reader {
		switch true {
		case *target && axtTargetGap(each):
			axt.WriteToFileHandle(writer, each, index)
			index++
		case *query && axtQueryGap(each):
			axt.WriteToFileHandle(writer, each, index)
			index++
		case *target && *query:
			log.Fatalf("Error: Must select gaps for target or query, not both...\n")
		}
	}
}

func axtTargetGap(record *axt.Axt) bool {
	log.SetFlags(0)
	if dna.CountBase(record.RSeq, dna.N) != 0 && dna.CountBase(record.QSeq, dna.N) == 0 {
		//log.Printf("%s", axt.ToString(record, 0))
		return true
	} else {
		return false
	}
}

func axtQueryGap(record *axt.Axt) bool {
	log.SetFlags(0)
	if dna.CountBase(record.QSeq, dna.N) != 0 && dna.CountBase(record.RSeq, dna.N) == 0 {
		//log.Printf("%s", axt.ToString(record, 0))
		return true
	} else {
		return false
	}
}

//TODO: Coming soon: axt alignment to Fasta

func axtToFa(input string, output string) {
	file := fileio.EasyOpen(input)
	reader := make(chan *axt.Axt)
	writer := fileio.EasyCreate(output)
	go axt.ReadToChan(file, reader)

	for each := range reader {
		fasta.WriteFasta(writer, axtSeq(each, nil), 50)
	}
}

//if target sequence contains Ns, uses query non N bases to fill Ns
func axtSeq(axtRecord *axt.Axt, fa *fasta.Fasta) *fasta.Fasta {
	concensus := &fasta.Fasta{
		Name: fmt.Sprintf("%s", axtRecord.RName),
		Seq:  make([]dna.Base, 0, len(fa.Seq)),
	}

	concensus.Seq = append(concensus.Seq, fa.Seq[:axtRecord.RStart-1]...)

	for i := 0; i < len(axtRecord.RSeq); i++ {
		if axtRecord.RSeq[i] == dna.N && axtRecord.QSeq[i] != dna.N {
			concensus.Seq = append(concensus.Seq, axtRecord.QSeq[i])
		} else {
			concensus.Seq = append(concensus.Seq, axtRecord.RSeq[i])
		}
	}

	concensus.Seq = append(concensus.Seq, fa.Seq[axtRecord.REnd:]...)

	if len(concensus.Seq) != len(fa.Seq) {
		log.Fatalf("Error: Sequence length is not the same...\n")
	}
	return concensus
}
