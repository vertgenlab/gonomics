// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faFindFast(inFile string, outFile string, windowSize int, chromName string, removeN bool, longWindowEfficient bool, verbose bool) {
	records := fasta.Read(inFile)

	if len(records) != 2 {
		log.Fatalf("Wrong number of sequences, expecting two, found %d.\n", len(records))
	}

	if len(records[0].Seq) != len(records[1].Seq) {
		log.Fatalf("Sequences are not of equal length")
	}

	if verbose {
		diff, denominator := countTotalDifference(records[0], records[1])
		log.Printf("I found %d total differences in %d aligned bases.\n", diff, denominator)
	}

	var bedList []bed.Bed
	if longWindowEfficient {
		file := fileio.EasyCreate(outFile)
		speedyWindowDifference(windowSize, records[0].Seq, records[1].Seq, chromName, removeN, verbose, file)
		err := file.Close()
		exception.PanicOnErr(err)
	} else {
		bedList = windowDifference(windowSize, records[0], records[1], chromName, removeN, verbose)
		bed.Write(outFile, bedList)
	}
}

func windowDifference(windowSize int, seq1 fasta.Fasta, seq2 fasta.Fasta, name string, removeN bool, verbose bool) []bed.Bed {
	var bedList []bed.Bed
	var referenceCounter = 0
	var reachedEnd bool = false
	var diff int = 0
	for alignmentCounter := 0; reachedEnd == false; alignmentCounter++ {
		if alignmentCounter%1000000 == 0 && verbose {
			log.Printf("alignmentCounter: %v\n", alignmentCounter)
		}

		if seq1.Seq[alignmentCounter] != dna.Gap {
			diff, reachedEnd = countWindowDifference(seq1, seq2, alignmentCounter, windowSize)
			if removeN && seqsContainN(seq1, seq2, alignmentCounter, windowSize) {
				referenceCounter++
			} else if !reachedEnd {
				current := bed.Bed{Chrom: name, ChromStart: referenceCounter,
					ChromEnd: referenceCounter + windowSize, Name: fmt.Sprintf("%d", referenceCounter), Score: diff, FieldsInitialized: 5}
				bedList = append(bedList, current)
				referenceCounter++
			}
		}
	}
	return bedList
}

func seqsContainN(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) bool {
	for i := start; i < start+windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == dna.N {
			return true
		} else if seq2.Seq[i] == dna.N {
			return true
		}
	}
	return false
}

func countWindowDifference(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) (int, bool) {
	diff := 0
	baseCount := 0
	var seq1Indel bool = false
	var seq2Indel bool = false
	var reachedEnd bool = false
	var i int = 0

	for i = start; baseCount < windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
				baseCount++
			}
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			baseCount++
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			baseCount++
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}

	if baseCount != windowSize {
		reachedEnd = true
	}
	return diff, reachedEnd
}

func countTotalDifference(seq1 fasta.Fasta, seq2 fasta.Fasta) (int, int) {
	var diff, denominator int
	var seq1Indel bool = false
	var seq2Indel bool = false

	for i := range seq1.Seq {
		if seq1.Seq[i] != dna.N || seq2.Seq[i] != dna.N {
			if seq1.Seq[i] == seq2.Seq[i] {
				seq1Indel = false
				seq2Indel = false
				denominator++
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
				denominator++
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
				denominator++
			} else {
				denominator++
			}
		}
	}
	return diff, denominator
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
	var chromName *string = flag.String("chrom", "", "Specify the chrom name")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var longWindowEfficient *bool = flag.Bool("longWindowEfficient", false, "Use the optimized algorithm for large window sizes.")
	var verbose *bool = flag.Bool("verbose", false, "Enable debug prints.")

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

	faFindFast(inFile, outFile, *windowSize, *chromName, *removeN, *longWindowEfficient, *verbose)
}
