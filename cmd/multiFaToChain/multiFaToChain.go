package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

type State byte

const (
	InAln  State = 0
	InTGap State = 1
	InQGap State = 2
)

func multiFaToChain(inFile string, tName string, qName string, outFile string, swapTandQ bool) {
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

	if swapTandQ {
		records[0], records[1] = records[1], records[0]
	}

	var answer []chain.Chain = make([]chain.Chain, 0)
	var header chain.HeaderComments = chain.HeaderComments{} //blank header to write
	var recordsCopy []fasta.Fasta = fasta.CopyAll(records)
	recordsCopy = fasta.RemoveGaps(recordsCopy)
	var currAlignment []chain.BaseStats = make([]chain.BaseStats, 0)
	var currBaseStats chain.BaseStats
	var prevState State
	var currState State
	var doubleGap bool

	//initialize state and first baseStats struct
	prevState, doubleGap = queryState(records, 0)
	currBaseStats = chain.BaseStats{Size: 0, TBases: 0, QBases: 0}

	for i := range records[0].Seq {
		currState, doubleGap = queryState(records, i)
		if doubleGap {
			continue
		}
		if prevState == currState { //if our state is the same at this position
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

	var TEnd = len(recordsCopy[0].Seq)
	var QEnd = len(recordsCopy[1].Seq)
	//if there are trailing Tbases or Qbases after the last alignment, we trim them away from the header
	//both T and Q are obligately positive stranded, so no need to check for strand.
	if currAlignment[len(currAlignment)-1].TBases > 0 {
		TEnd -= currAlignment[len(currAlignment)-1].TBases
	}

	if currAlignment[len(currAlignment)-1].QBases > 0 {
		QEnd -= currAlignment[len(currAlignment)-1].QBases
	}

	var currAnswer chain.Chain = chain.Chain{
		Score:     100, //dummy value
		TName:     tName,
		TSize:     len(recordsCopy[0].Seq),
		TStrand:   true,
		TStart:    0,
		TEnd:      TEnd,
		QName:     qName,
		QSize:     len(recordsCopy[1].Seq),
		QStrand:   true,
		QStart:    0,
		QEnd:      QEnd,
		Alignment: currAlignment,
		Id:        1,
	}

	answer = append(answer, currAnswer)
	chain.Write(outFile, answer, header)
}

// queryState determines the State of a multiFa alignment at a particular input index position. Second bool returns true
// if a doubleGap is observed.
func queryState(records []fasta.Fasta, index int) (State, bool) {
	if dna.DefineBase(records[0].Seq[index]) || records[0].Seq[index] == dna.N || records[0].Seq[index] == dna.LowerN {
		if dna.DefineBase(records[1].Seq[index]) || records[1].Seq[index] == dna.N || records[1].Seq[index] == dna.LowerN {
			return InAln, false
		} else if records[1].Seq[index] == dna.Gap {
			return InQGap, false
		} else {
			log.Fatalf("Unrecognized dna base in the query sequence: %s.", dna.BaseToString(records[1].Seq[index]))
		}
	} else if records[0].Seq[index] == dna.Gap {
		if dna.DefineBase(records[1].Seq[index]) || records[1].Seq[index] == dna.N || records[1].Seq[index] == dna.LowerN {
			return InTGap, false
		} else if records[1].Seq[index] == dna.Gap {
			return InAln, true //inAln is a dummy here, we care about the bool return signaling a doubleGap in the alignment.
		}
	} else {
		log.Fatalf("Unrecognized dna base in the reference sequence: %s.", dna.BaseToString(records[0].Seq[index]))
		return 3, true
	}
	return 3, true
}

func usage() {
	fmt.Print(
		"multiFaToChain - Convert a pairwise multiFa format alignment to a chain file. First species is the target by default.\n" +
			"Usage:\n" +
			" multiFaToChain input.fa tName qName output.chain\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var swapTandQ *bool = flag.Bool("swapTandQ", false, "Swap the target and query in the output chain file.")

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

	multiFaToChain(inFile, tName, qName, outFile, *swapTandQ)
}
