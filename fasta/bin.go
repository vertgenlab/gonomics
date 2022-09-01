package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
)

func BinFasta(genome []Fasta, binNum int) map[int][]Fasta {
	var remainingBases, totalBases, lastPos int
	var lastName string
	var answer map[int][]Fasta
	for c := range genome {
		totalBases = totalBases + len(genome[c].Seq)
	}
	//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with baseCap,
	//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
	//	could do a check to see if remainder is super large and then handle as necessary
	equalBinNum := binNum - 1
	baseCap := totalBases / equalBinNum

	for j := range answer {
		for i := range genome {
			var done int
			var currFa Fasta
			var currBases []dna.Base
			if j == binNum { //we have finished making all bins but the last bin which will hold the remainder of bases
				answer[j], done = fillLastBin(genome, i)
				i = done
			}

			value, exists := answer[j]
			if value != nil { //old bin
				filledBases := calcNumBasesInBin(answer[j])
				if filledBases < baseCap { //bin not full
					remainingBases = baseCap - filledBases
					currBases = appendRangeOfBases(genome[i].Seq, 0, remainingBases-1)
				}
			} else if !exists || value == nil { //new bin
				lastName, lastPos = findLastBinnedPos(answer[j-1])
				if lastName == genome[i].Name && lastPos < len(genome[i].Seq) { //new bin and old contig
					remainingBases = len(genome[i].Seq) - lastPos
					currBases = appendRangeOfBases(genome[i].Seq, lastPos+1, lastPos+remainingBases)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
				} else if lastName == genome[i].Name && lastPos == len(genome[i].Seq) { //new bin and new contig
					currBases = appendRangeOfBases(genome[i].Seq, 0, baseCap)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
				}
			}
		}
	}
	return answer
}

func calcNumBasesInBin(f []Fasta) int {
	totalBases := 0
	for i := range f {
		totalBases = totalBases + len(f[i].Seq)
	}
	return totalBases
}

func findLastBinnedPos(f []Fasta) (string, int) {
	var lastPos int
	var lastFasta string

	for i := range f {
		lastFasta = f[i].Name
		lastPos = len(f[i].Seq)
	}

	return lastFasta, lastPos
}

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
			if b >= start && b <= stop {
				answer = append(answer, bases[b])
			}
		}
	}

	return answer
}

func fillLastBin(f []Fasta, currRec int) (ans []Fasta, final int) {
	var answer []Fasta
	var bases []dna.Base
	var i int

	for i = currRec; i < len(f); i++ {
		bases = appendRangeOfBases(bases, 0, len(f[i].Seq))
	}

	return answer, i
}
