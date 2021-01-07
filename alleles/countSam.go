package alleles

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

//TODO: Merge with countGiraf functions using interfaces???
// GoCountSamAlleles is a wrapper for CountSamAlleles that manages channel closure.
func GoCountSamAlleles(samFilename string, reference []*fasta.Fasta, minMapQ int64) <-chan *Allele {
	answer := make(chan *Allele)
	var wg sync.WaitGroup
	wg.Add(1)
	go CountSamAlleles(answer, samFilename, reference, minMapQ, &wg)

	go func() {
		wg.Wait()
		close(answer)
	}()

	return answer
}

// CountSamAlleles counts the alleles in a sam file aligned to a linear reference (fasta) and sends them to an input channel sending the allele count for each position in the reference covered by the sam file
func CountSamAlleles(answer chan<- *Allele, samFilename string, reference []*fasta.Fasta, minMapQ int64, wg *sync.WaitGroup) {
	samChan, _ := sam.GoReadToChan(samFilename)
	var currAlleles = make(map[Coordinate]*AlleleCount)
	var runningCount = make([]*Coordinate, 0)
	var progress int // TODO: Make option to print progress

	fasta.AllToUpper(reference)
	ref := fasta.FastaMap(reference)

	for read := range samChan {
		runningCount = sendPassedPositionsSam(answer, read, samFilename, runningCount, currAlleles)
		currAlleles, runningCount = countSamRead(read, currAlleles, runningCount, ref, minMapQ, progress)
	}
	wg.Done()
}

// sendPassedPositionsSam sends positions that have been passed in the file
func sendPassedPositionsSam(answer chan<- *Allele, aln *sam.SamAln, samFilename string, runningCount []*Coordinate, currAlleles map[Coordinate]*AlleleCount) []*Coordinate {
	for i := 0; i < len(runningCount); i++ {

		if runningCount[i].Chr != aln.RName {
			answer <- &Allele{samFilename, currAlleles[*runningCount[i]], runningCount[i]}
			delete(currAlleles, *runningCount[i])

			// Catch instance where every entry in running count is sent
			// Delete all of runningCount
			if i == len(runningCount)-1 {
				runningCount = nil
			}
			continue
		}

		if runningCount[i].Pos < int(aln.Pos-1) {
			answer <- &Allele{samFilename, currAlleles[*runningCount[i]], runningCount[i]}
			delete(currAlleles, *runningCount[i])

			// Catch instance where every entry in running count is sent
			// Delete all of runningCount
			if i == len(runningCount)-1 {
				runningCount = nil
			}

		} else {
			// Remove sent values from count
			runningCount = runningCount[i:]
			break
		}
	}
	return runningCount
}

// countSamRead adds the bases in a single sam read to the currAlleles map
func countSamRead(aln *sam.SamAln, currAlleles map[Coordinate]*AlleleCount, runningCount []*Coordinate, ref map[string][]dna.Base, minMapQ int64, progress int) (map[Coordinate]*AlleleCount, []*Coordinate) {

	if aln.Cigar[0].Op == '*' {
		return currAlleles, runningCount
	}

	if aln.MapQ < minMapQ {
		return currAlleles, runningCount
	}

	var RefIndex, SeqIndex, OrigRefIndex int
	var currentSeq []dna.Base
	var i, j, k int
	var currentIndel Indel
	var indelSeq []dna.Base
	var Match bool

	// Count the bases
	progress++
	SeqIndex = 0
	RefIndex = int(aln.Pos - 1)

	for i = 0; i < len(aln.Cigar); i++ {
		currentSeq = aln.Seq

		if aln.Cigar[i].Op == 'D' {
			OrigRefIndex = RefIndex
			indelSeq = make([]dna.Base, 1)

			// First base in indel is the base prior to the indel sequence per VCF standard format
			indelSeq[0] = ref[aln.RName][OrigRefIndex-1]

			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {

				// If the position has already been added to the map, move along
				_, ok := currAlleles[Coordinate{aln.RName, RefIndex}]

				// If the position is NOT in the map, initialize
				if !ok {
					currAlleles[Coordinate{aln.RName, RefIndex}] = &AlleleCount{
						Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
					runningCount = append(runningCount, &Coordinate{aln.RName, RefIndex})
				}

				// Keep track of deleted sequence
				indelSeq = append(indelSeq, ref[aln.RName][RefIndex])

				currAlleles[Coordinate{aln.RName, RefIndex}].Counts++
				RefIndex++
			}

			Match = false
			for j = 0; j < len(currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel); j++ {
				// If the deletion has already been seen before, increment the existing entry
				// For a deletion the indelSeq should match the Ref
				if dna.CompareSeqsIgnoreCase(indelSeq, currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel[j].Ref) == 0 &&
					dna.CompareSeqsIgnoreCase(indelSeq[:1], currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel[j].Alt) == 0 {
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel[j].CountF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel[j].CountR++
					}

					Match = true
					break
				}
			}

			// If the deletion has not been seen before, then append it to the Del slice
			// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
			if Match == false {

				currentIndel = Indel{indelSeq, indelSeq[:1], 0, 0}
				if sam.IsForwardRead(aln) == true {
					currentIndel.CountF++
				} else if sam.IsReverseRead(aln) == false {
					currentIndel.CountR++
				}
				currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel = append(currAlleles[Coordinate{aln.RName, OrigRefIndex}].Indel, currentIndel)
			}

			//Handle insertion relative to ref
			//The base after the inserted sequence is annotated with an Ins read
		} else if aln.Cigar[i].Op == 'I' {

			// If the position has already been added to the map, move along
			_, ok := currAlleles[Coordinate{aln.RName, RefIndex}]

			// If the position is NOT in the map, initialize
			if !ok {
				currAlleles[Coordinate{aln.RName, RefIndex}] = &AlleleCount{
					Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
				runningCount = append(runningCount, &Coordinate{aln.RName, RefIndex})
			}

			// Loop through read sequence and keep track of the inserted bases
			indelSeq = make([]dna.Base, 1)

			// First base in indel is the base prior to the indel sequence per VCF standard format
			indelSeq[0] = ref[aln.RName][RefIndex-1]

			// Keep track of inserted sequence by moving along the read
			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {
				indelSeq = append(indelSeq, currentSeq[SeqIndex])
				SeqIndex++
			}

			Match = false
			for j = 0; j < len(currAlleles[Coordinate{aln.RName, RefIndex}].Indel); j++ {
				// If the inserted sequence matches a previously inserted sequence, then increment the count
				// For an insertion, the indelSeq should match the Alt
				if dna.CompareSeqsIgnoreCase(indelSeq, currAlleles[Coordinate{aln.RName, RefIndex}].Indel[j].Alt) == 0 &&
					dna.CompareSeqsIgnoreCase(indelSeq[:1], currAlleles[Coordinate{aln.RName, RefIndex}].Indel[j].Ref) == 0 {
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].Indel[j].CountF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].Indel[j].CountR++
					}
					Match = true
					break
				}
			}

			if Match == false {
				currentIndel = Indel{indelSeq[:1], indelSeq, 0, 0}
				if sam.IsForwardRead(aln) == true {
					currentIndel.CountF++
				} else if sam.IsReverseRead(aln) == true {
					currentIndel.CountR++
				}
				currAlleles[Coordinate{aln.RName, RefIndex}].Indel = append(currAlleles[Coordinate{aln.RName, RefIndex}].Indel, currentIndel)
			}

			// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

			//Handle matching pos relative to ref
		} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {

				//if the position has already been added to the matrix, move along
				_, ok := currAlleles[Coordinate{aln.RName, RefIndex}]

				//if the position is NOT in the matrix, add it
				if !ok {
					currAlleles[Coordinate{aln.RName, RefIndex}] = &AlleleCount{
						Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
					runningCount = append(runningCount, &Coordinate{aln.RName, RefIndex})
				}

				switch currentSeq[SeqIndex] {
				case dna.A:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseAF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseAR++
					}
					currAlleles[Coordinate{aln.RName, RefIndex}].Counts++
				case dna.T:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseTF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseTR++
					}
					currAlleles[Coordinate{aln.RName, RefIndex}].Counts++
				case dna.G:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseGF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseGR++
					}
					currAlleles[Coordinate{aln.RName, RefIndex}].Counts++
				case dna.C:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseCF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Coordinate{aln.RName, RefIndex}].BaseCR++
					}
					currAlleles[Coordinate{aln.RName, RefIndex}].Counts++
				}
				SeqIndex++
				RefIndex++
			}
		} else if aln.Cigar[i].Op != 'H' {
			SeqIndex = SeqIndex + int(aln.Cigar[i].RunLength)
		}
	}
	return currAlleles, runningCount
}
