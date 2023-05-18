package fasta

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
)

// BinGenomeNoBreaks takes in an entire genome which is sorted largest to smallest contig and breaks up the fasta so
// that smaller contigs get combined into a single fasta, while large contigs become a single fasta on their own. The user
// must specify the number of bins for the genome to be broken into, the genome must have more contigs than bins in order
// to combine any contigs and equal number bins to contigs if each contig gets its own record. The bins will all be
// filled with the first contig encountered when it's empty, and then the smallest of those bins will be filled when the
// contig is equal to binNum+1. The minSize option allows for a user to specify a minimum length of sequence to go into
// each bin and in this case the number of bins returned depends on the minSize and the binNum will be ignored.
func BinGenomeNoBreaks(genome []Fasta, binNum int, minSize int) map[int][]Fasta {
	var bins map[int][]Fasta

	if minSize != -1 {
		bins = binMinSize(genome, minSize)
	} else {
		bins = make(map[int][]Fasta, binNum)
		if len(genome) > binNum {
			for n := 0; n < binNum; n++ {
				bins[n] = append(bins[n], genome[n])
			}
			bins = fillSmallestBin(bins, genome)
		} else if len(genome) == binNum {
			for n := 0; n < binNum; n++ {
				bins[n] = append(bins[n], genome[n])
			}
		} else {
			log.Fatal("Number of bins is greater than the number of contigs in the given genome. Reduce bin number.")
		}
	}
	return bins
}

// fillSmallestBin fills the bin found by findSmallestBin and adds the next contig to it.
func fillSmallestBin(bins map[int][]Fasta, genome []Fasta) map[int][]Fasta {
	for i := len(bins); i < len(genome); i++ {
		b := findSmallestBin(bins)
		bins[b] = append(bins[b], genome[i])
	}

	return bins
}

// findSmallestBin finds the bin with the least sequence and returns which bin it is so that it can be filled with the next contig.
func findSmallestBin(bins map[int][]Fasta) int {
	var smallest int
	var sizeSmallest int
	for i := 0; i < len(bins); i++ {
		rec, exists := bins[i]
		if !exists {
			log.Panic("Map was not filled properly in BinGenomeNoBreaks.")
		}
		size := calcNumBasesInBin(rec)
		if size < sizeSmallest {
			sizeSmallest = size
			smallest = i
		} else if sizeSmallest == 0 {
			sizeSmallest = size
		}
	}
	return smallest
}

// binMinSize fills bins to a minimum length of sequence as specified by the user and returns whatever number of bins it may make.
func binMinSize(genome []Fasta, min int) map[int][]Fasta {
	var length int
	bins := make(map[int][]Fasta, len(genome))

	for i := range genome {
		if len(bins) == 0 {
			bins[0] = append(bins[0], genome[i])
		} else if len(bins) > 0 {
			if len(genome[i].Seq) > min {
				length = len(bins)
				for j := len(bins); j < length+1; j++ {
					value, ok := bins[j]
					if !ok || value == nil {
						bins[j] = append(bins[j], genome[i])
					}
				}
			} else {
				k := findBinBelowMin(bins, min)
				if k < 0 && i+1 == len(genome) {
					bins[len(bins)-1] = append(bins[len(bins)-1], genome[i])
				} else if k < 0 {
					value := bins[len(bins)]
					if value == nil {
						bins[len(bins)] = append(bins[len(bins)], genome[i])
					}
				} else {
					bins[k] = append(bins[k], genome[i])
				}
			}
		}
	}

	return bins
}

// findBinBelowMin finds the first bin that is filled below the minimum length of sequence and returns which bin it is
// so that it can be filled with the next contig below the minimum length.
func findBinBelowMin(bins map[int][]Fasta, min int) int {
	answer := -1

	for i := 0; i < len(bins); i++ {
		bases := calcNumBasesInBin(bins[i])
		if bases < min {
			answer = i
		}
	}
	return answer
}

// calcNumBasesInBin will determine the number of bases that already exist in any given bin.
func calcNumBasesInBin(f []Fasta) int {
	totalBases := 0
	for i := range f {
		totalBases = totalBases + len(f[i].Seq)
	}
	return totalBases
}

// BinFasta takes in a slice of fastas and breaks it up into x number of fastas with relatively
// equal sequence in each, where x equals the number of bins specified.
func BinFasta(genome []Fasta, binNum int) map[int][]Fasta {
	if binNum == 0 {
		log.Panic("Number of bins must be greater than zero.")
	}

	var totalBases, lastPos, baseCap, remainder, splitBases, equalBinNum int
	var lastName string
	answer := make(map[int][]Fasta, binNum)
	for c := range genome {
		totalBases = totalBases + len(genome[c].Seq)
	}

	equalBinNum = binNum - 1
	remainder = totalBases % binNum
	if remainder == 0 {
		baseCap = totalBases / binNum
	} else {
		splitBases = totalBases - remainder //number of bases to split between equalBins
		baseCap = splitBases / equalBinNum
	}

	for j := 0; j < binNum; j++ {
		for i := 0; i < len(genome); i++ {
			var remainingBases int
			var done int
			var currFa Fasta
			var currBases []dna.Base

			if j == binNum { //we have finished making all bins but the last bin which will hold the remainder of bases
				answer[j], done = fillLastBin(genome, i)
				i = done
			}

			value, exists := answer[j]
			if j == 0 && i == 0 {
				currBases = appendRangeOfBases(genome[i].Seq, 0, baseCap)
				currFa.Name = genome[i].Name
				currFa.Seq = currBases
				answer[j] = append(answer[j], currFa)
				lastPos = baseCap
				continue
			} else if value != nil { //old bin
				filledBases := calcNumBasesInBin(answer[j])
				if filledBases < baseCap { //bin not full
					remainingBases = baseCap - filledBases
					currBases = appendRangeOfBases(genome[i].Seq, 0, remainingBases)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
					lastPos = len(currBases)
				} else {
					continue
				}
			} else if !exists || value == nil { //new bin
				lastName = findLastBinnedName(answer[j-1])
				if lastName == genome[i].Name && lastPos < len(genome[i].Seq) { //new bin and old contig
					remainingBases = len(genome[i].Seq) - lastPos
					currBases = appendRangeOfBases(genome[i].Seq, lastPos, lastPos+remainingBases)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
					lastPos = lastPos + len(currBases)
				} else if lastName == genome[i].Name && lastPos == len(genome[i].Seq) { //new bin and new contig
					lastPos = 0
				} else if lastName != "" && lastPos == 0 {
					currBases = appendRangeOfBases(genome[i].Seq, 0, baseCap)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
					lastPos = baseCap
				}
			}
		}
	}
	return answer
}

// findLastBinnedName determines the last contig that was handled in the given bin (f).
func findLastBinnedName(f []Fasta) string {
	var lastFasta string

	for i := range f {
		lastFasta = f[i].Name
	}

	return lastFasta
}

// appendRangeOfBases will create a slice of bases ranging from the start to stop positions in the records given (bases).
func appendRangeOfBases(bases []dna.Base, start int, stop int) []dna.Base {
	var answer []dna.Base

	if len(bases) < stop {
		for b := range bases {
			if b >= start && b <= len(bases) {
				answer = append(answer, bases[b])
			}
		}
	} else {
		for b := range bases {
			if b >= start && b <= stop-1 {
				answer = append(answer, bases[b])
			}
		}
	}

	return answer
}

// fillLastBin takes all remaining contigs once BinFasta has reached the last bin and puts them all in the final bin regardless of sequence length.
func fillLastBin(f []Fasta, currRec int) (ans []Fasta, final int) {
	var answer []Fasta
	var bases []dna.Base
	var i int

	for i = currRec; i < len(f); i++ {
		bases = appendRangeOfBases(bases, 0, len(f[i].Seq))
	}

	return answer, i
}
