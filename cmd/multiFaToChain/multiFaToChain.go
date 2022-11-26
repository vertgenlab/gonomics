package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

type State byte

const (
	InAln  State = 0
	InTGap State = 1
	InQGap State = 2
)

func multiFaToChain(inFile string, tName string, qName string, outFile string) {
	var records []fasta.Fasta = fasta.Read(inFile)

	//preflight checks
	if len(records) != 2 {
		log.Fatalf("multiFaToChain accepts only pairwise multiFa alignment files.")
	}
	if len(records[0].Seq) != len(records[1].Seq) {
		log.Fatalf("Both sequences must be of the same alignment length.")
	}
	if len(records[0].Seq) < 1 {
		log.Fatalf("MultiFaToChain expects non-empty DNA sequences.")
	}

	var answer []chain.Chain = make([]chain.Chain, 0)
	var header chain.HeaderComments = chain.HeaderComments{} //blank header to write
	var recordsCopy []fasta.Fasta = fasta.CopyAll(records)
	recordsCopy = fasta.RemoveGaps(recordsCopy)
	var currAlignment []chain.BaseStats = make([]chain.BaseStats, 0)
	var currBaseStats chain.BaseStats
	var prevState State
	var currState State

	//initialize state and first baseStats struct
	prevState = queryState(records, 0, 3) //previous state hardcoded to 3.
	currBaseStats = chain.BaseStats{Size: 0, TBases: 0, QBases: 0}

	for i := range records[0].Seq {
		currState = queryState(records, i, prevState)
		if prevState == queryState(records, i, prevState) { //if our state is the same at this position
			switch prevState {
			case InAln:
				currBaseStats.Size++
			case InTGap:
				currBaseStats.QBases++
			case InQGap:
				currBaseStats.TBases++
			}
		} else { //our state has changed, so we can might have to write out the current line
			switch prevState {
			case InAln: //we leave an alignment and enter either a Tgap or Qgap
				if currState == InQGap {
					currBaseStats.TBases++
					prevState = InQGap
				} else if currState == InTGap {
					currBaseStats.QBases++
					prevState = InTGap
				} else {
					log.Fatalf("Unrecognized state.")
				}
			case InTGap:
				if currState == InAln { //if we go from a T gap to alignment, we write out the old baseStats and start a new block.
					currAlignment = append(currAlignment, currBaseStats)
					currBaseStats = chain.BaseStats{Size: 1, TBases: 0, QBases: 0} //size is 1 to count current position
					prevState = InAln
				} else if currState == InQGap { //I'm not positive this is allowed, but if we have a t gap that switches to a q gap we'll make it part of the same line.
					prevState = InQGap
					currBaseStats.TBases++
				} else {
					log.Fatalf("Unrecognized state.")
				}
			case InQGap:
				if currState == InAln {
					currAlignment = append(currAlignment, currBaseStats)
					currBaseStats = chain.BaseStats{Size: 1, TBases: 0, QBases: 0}
					prevState = InAln
				} else if currState == InTGap { //again, not certain this can be allowed.
					prevState = InTGap
					currBaseStats.QBases++
				}
			default:
				log.Fatalf("Unrecognized state.")
			}
		}
	}
	currAlignment = append(currAlignment, currBaseStats)

	var currAnswer chain.Chain = chain.Chain{
		Score:     100, //dummy value
		TName:     tName,
		TSize:     len(recordsCopy[0].Seq),
		TStrand:   true,
		TStart:    0,
		TEnd:      len(recordsCopy[0].Seq),
		QName:     qName,
		QSize:     len(recordsCopy[1].Seq),
		QStrand:   true,
		QStart:    0,
		QEnd:      len(recordsCopy[1].Seq),
		Alignment: currAlignment,
		Id:        1,
	}

	answer = append(answer, currAnswer)
	chain.Write(outFile, answer, header)
}

func queryState(records []fasta.Fasta, index int, prevState State) State {
	if dna.DefineBase(records[0].Seq[index]) || records[0].Seq[index] == dna.N || records[0].Seq[index] == dna.LowerN {
		if dna.DefineBase(records[1].Seq[index]) || records[1].Seq[index] == dna.N || records[1].Seq[index] == dna.LowerN {
			return InAln
		} else if records[1].Seq[index] == dna.Gap {
			return InQGap
		} else {
			log.Fatalf("Unrecognized dna base in the query sequence: %s.", dna.BaseToString(records[1].Seq[index]))
		}
	} else if records[0].Seq[index] == dna.Gap {
		if dna.DefineBase(records[1].Seq[index]) || records[1].Seq[index] == dna.N || records[1].Seq[index] == dna.LowerN {
			return InTGap
		} else if records[1].Seq[index] == dna.Gap {
			return prevState
		}
	} else {
		log.Fatalf("Unrecognized dna base in the reference sequence: %s.", dna.BaseToString(records[0].Seq[index]))
		return 3
	}
	return 3
}

func usage() {
	fmt.Print(
		"multiFaToChain - Convert a pairwise multiFa format alignment to a chain file.\n" +
			"Usage:\n" +
			" multiFaToChain input.fa tName qName output.chain\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	tName := flag.Arg(1)
	qName := flag.Arg(2)
	outFile := flag.Arg(3)

	multiFaToChain(inFile, tName, qName, outFile)
}
