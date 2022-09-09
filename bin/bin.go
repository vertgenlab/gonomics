package bin

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

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
