package bin

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

//BinGenomeNoBreaks takes in an entire genome which is sorted largest to smallest contig and breaks up the fasta so
//that smaller contigs get combined into a single fasta, while large contigs become a single fasta on their own. The user
//must specify the number of bins for the genome to be broken into, the genome must have more contigs than bins in order
//to combine any contigs and equal number bins to contigs if each contig gets its own record. The bins will all be
//filled with the first contig encountered when it's empty, and then the smallest of those bins will be filled when the
//contig is equal to binNum+1. The minSize option allows for a user to specify a minimum length of sequence to go into
//each bin and in this case the number of bins returned depends on the minSize and the binNum will be ignored.
func BinGenomeNoBreaks(genome []fasta.Fasta, binNum int, minSize int) map[int][]fasta.Fasta {
	bins := make(map[int][]fasta.Fasta, binNum)
	if minSize != -1 {
		bins = binMinSize(genome, minSize)
	} else {
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

//fillSmallestBin fills the bin found by findSmallestBin and adds the next contig to it.
func fillSmallestBin(bins map[int][]fasta.Fasta, genome []fasta.Fasta) map[int][]fasta.Fasta {

	for i := len(bins); i < len(genome); i++ {
		b := findSmallestBin(bins)
		bins[b] = append(bins[b], genome[i])
	}

	return bins
}

//findSmallestBin finds the bin with the least sequence and returns which bin it is so that it can be filled with the next contig.
func findSmallestBin(bins map[int][]fasta.Fasta) int {
	var smallest int
	var sizeSmallest int
	for i := range bins {
		rec, exists := bins[i]
		if !exists {
			log.Panic("Map was not filled properly in BinGenomeNoBreaks.")
		}
		size := calcNumBasesInBin(rec)
		if size < sizeSmallest {
			sizeSmallest = size
			smallest = i
		}
	}
	return smallest
}

//binMinSize fills bins to a minimum length of sequence as specified by the user and returns whatever number of bins it may make.
func binMinSize(genome []fasta.Fasta, min int) map[int][]fasta.Fasta {
	var bins map[int][]fasta.Fasta

	for i, chr := range genome {
		if len(chr.Seq) > min {
			for j := range bins {
				value, ok := bins[i]
				if !ok || value == nil {
					bins[j] = append(bins[j], chr)
				}
			}
		} else {
			k := findBinBelowMin(bins, min)
			if k < 0 {
				for j := range bins {
					value, ok := bins[i]
					if !ok || value == nil {
						bins[j] = append(bins[j], chr)
					}
				}
			} else {
				bins[k] = append(bins[k], chr)
			}
		}
	}

	return bins
}

//findBinBelowMin finds the first bin that is filled below the minimum length of sequence and returns which bin it is
//so that it can be filled with the next contig below the minimum length
func findBinBelowMin(bins map[int][]fasta.Fasta, min int) int {
	answer := -1

	for i := range bins {
		fast, _ := bins[i]
		for f := range fast {
			if len(fast[f].Seq) < min {
				answer = i
				return answer
			}
		}
	}
	return answer
}

//TODO: change all logic below here to break a single chrom into whatever number of bins

//BinFasta takes in a slice of fastas and breaks it up into x number of fastas with relatively
//equal sequence in each, where x equals the number of bins specified
func BinFasta(genome []fasta.Fasta, binNum int) map[int][]fasta.Fasta {
	if binNum == 0 {
		log.Panic("Number of bins must be greater than zero.")
	}

	var totalBases, lastPos, baseCap, remainder, splitBases, equalBinNum int
	var lastName string
	answer := make(map[int][]fasta.Fasta, binNum)
	for c := range genome {
		totalBases = totalBases + len(genome[c].Seq)
	}

	//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with baseCap,
	//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
	//	could do a check to see if remainder is super large and then handle as necessary
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
			var currFa fasta.Fasta
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

//calcNumBasesInBin will determine the number of bases that already exist in any given bin
func calcNumBasesInBin(f []fasta.Fasta) int {
	totalBases := 0
	for i := range f {
		totalBases = totalBases + len(f[i].Seq)
	}
	return totalBases
}

//findLastBinnedName determines the last contig that was handled in the given bin (f)
func findLastBinnedName(f []fasta.Fasta) string {
	var lastFasta string

	for i := range f {
		lastFasta = f[i].Name
	}

	return lastFasta
}

//appendRangeOfBases will create a slice of bases ranging from the start to stop positions in the records given (bases)
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

//fillLastBin takes all remaining contigs once BinFasta has reached the last bin and puts them all in the final bin regardless of sequence length
func fillLastBin(f []fasta.Fasta, currRec int) (ans []fasta.Fasta, final int) {
	var answer []fasta.Fasta
	var bases []dna.Base
	var i int

	for i = currRec; i < len(f); i++ {
		bases = appendRangeOfBases(bases, 0, len(f[i].Seq))
	}

	return answer, i
}
