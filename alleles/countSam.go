package alleles

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

// GoSamPileup reads a sam file and streams per-base pileup info to the output channel.
func GoSamPileup(samFilename string, reference []fasta.Fasta, minMapQ uint8) <-chan Allele {
	answer := make(chan Allele, 1000)
	go SamPileup(answer, samFilename, reference, minMapQ)
	return answer
}

// SamPileup counts the alleles in a sam file aligned to a linear reference (fasta) and sends
// them to an input channel sending the allele count for each position in the reference covered by the sam file
func SamPileup(data chan<- Allele, samFilename string, reference []fasta.Fasta, minMapQ uint8) {
	samChan, _ := sam.GoReadToChan(samFilename)
	var currAlleles = make(Pileup)
	var currOpenCoords = make([]Coordinate, 0)

	ref := fasta.ToMap(reference)

	for read := range samChan {
		currOpenCoords = sendPassedPositionsSam(data, read, samFilename, currOpenCoords, currAlleles)
		currAlleles, currOpenCoords = countSamRead(read, currAlleles, currOpenCoords, ref, minMapQ)
	}
	close(data)
}

// sendPassedPositionsSam sends positions that have been passed in the file
func sendPassedPositionsSam(answer chan<- Allele, aln sam.Sam, currOpenCoords []Coordinate, currAlleles Pileup) []Coordinate {
	for i := 0; i < len(currOpenCoords); i++ {

		if currOpenCoords[i].Chr != aln.RName {
			answer <- &Allele{samFilename, currAlleles[*currOpenCoords[i]], currOpenCoords[i]}
			delete(currAlleles, *currOpenCoords[i])

			// Catch instance where every entry in running count is sent
			// Delete all of currOpenCoords
			if i == len(currOpenCoords)-1 {
				currOpenCoords = nil
			}
			continue
		}

		if currOpenCoords[i].Pos < int(aln.Pos-1) {
			answer <- &Allele{samFilename, currAlleles[*currOpenCoords[i]], currOpenCoords[i]}
			delete(currAlleles, *currOpenCoords[i])

			// Catch instance where every entry in running count is sent
			// Delete all of currOpenCoords
			if i == len(currOpenCoords)-1 {
				currOpenCoords = nil
			}

		} else {
			// Remove sent values from count
			currOpenCoords = currOpenCoords[i:]
			break
		}
	}
	return currOpenCoords
}

// countSamRead adds the bases in a single sam to the pileup.
func countSamRead(aln sam.Sam, currAlleles Pileup, currOpenCoords []Coordinate, ref fasta.FastaMap, minMapQ uint8) (Pileup, []Coordinate) {

	if aln.Cigar[0].Op == '*' {
		return currAlleles, currOpenCoords
	}

	if aln.MapQ < minMapQ {
		return currAlleles, currOpenCoords
	}

	var seqPos int
	currCoord := Coordinate{Chr: aln.RName, Pos: int(aln.Pos) - 1} // sam is native 1-base

	for _, cig := range aln.Cigar {
		switch cig.Op {
		case 'M':
		case 'I':
		case 'D':
		case 'N':
		case 'S':
		case 'H':
		case 'P':
		case 'X':
		case '=':
		}
	}

	return currAlleles, currOpenCoords
}
