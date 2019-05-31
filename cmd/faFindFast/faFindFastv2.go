package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faFindFast(inFile string, outFile string, windowSize int) {
	records, err := fasta.Read(inFile)
	if err != nil {
		log.Fatal(err)
	}

	if len(records) != 2 {
		log.Fatalf("Wrong number of sequences, expecting two, found %d.\n", len(records))
	}

	if len(records[0].Seq) != len(records[1].Seq) {
		log.Fatalf("Sequences are not of equal length")
	}

	diff := countTotalDifference(records[0], records[1])
	fmt.Printf("I found %d total differences.\n", diff)

	bedList := windowDifference(windowSize, records[0], records[1])
	bed.Write(outFile, bedList, 5)
}

func windowDifference(windowSize int, seq1 fasta.Fasta, seq2 fasta.Fasta) []*bed.Bed {

	var bedList []*bed.Bed

	for i := 0; i < (len(seq1.Seq) - (windowSize - 1)); i++ {
		current := bed.Bed{Chrom: seq1.Name, ChromStart: int64(i), ChromEnd: int64(i + windowSize), Name: fmt.Sprintf("%d", i), Score: 0}
		current.Score = int64(countWindowDifference(seq1, seq2, i, windowSize))
		bedList = append(bedList, &current)
	}

	return bedList
}

func countWindowDifference(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) int {
	diff := 0
	var seq1Indel bool = false
	var seq2Indel bool = false

	for i := start; i < (start + windowSize); i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
			}
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}
	return diff
}

func countTotalDifference(seq1 fasta.Fasta, seq2 fasta.Fasta) int {
	diff := 0
	var seq1Indel bool = false
	var seq2Indel bool = false

	for i, _ := range seq1.Seq {
		if seq1.Seq[i] == seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			continue
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false

			if seq1Indel {
				continue
			} else {
				seq1Indel = true
				diff++
				continue
			}
		} else if seq2.Seq[i] == dna.Gap {
			seq1Indel = false

			if seq2Indel {
				continue
			} else {
				seq2Indel = true
				diff++
				continue
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			diff++
		}
	}
	return diff
}

func usage() {
	fmt.Print(
		"faFindFast - Returns regions with highest SNP density\n" +
			"Usage:\n" +
			" faFindFast input.fa output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
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

	faFindFast(inFile, outFile, *windowSize)
}
