package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
)

func BinFasta(fa []Fasta, binNum int) {
	var totalBases int
	bins := make(map[int][]Fasta, binNum)
	contigs := len(fa)

	if contigs < binNum {
		answer := fewerContigsThanBins(fa, binNum)
		totalBases = 0
		for f := range fa {
			totalBases = totalBases + len(fa[f].Seq)
		}
		//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with bases per equal bin,
		//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
		//	could do a check to see if remainder is super large and then handle as necessary
		equalBinNum := binNum - 1
		basesPerEqualBin := totalBases / equalBinNum
		leftOvers := totalBases % equalBinNum

	} else if contigs == binNum {
		//TODO: split first 10-15 in half and add to last 10-15
	} else if contigs > binNum {

	}

}

func fewerContigsThanBins(genome []Fasta, binNum int) map[int][]Fasta {
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
	leftOver := totalBases % equalBinNum

	for j := range answer {
		//TODO: have a check for when we are finished with all equal bins
		for i := range genome {
			var currFa Fasta
			var currBases []dna.Base
			value, exists := answer[j]
			if value != nil { //old bin
				var filledBases int
				for s := range answer[j] {
					filledBases = filledBases + len(answer[j][s].Seq)
				}
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

	//for k := range answer {
	//	existingBases := 0
	//	if answer[k] != nil {
	//		existingBases = calcNumBasesInBin(answer[k])
	//		if existingBases < baseCap {
	//			remainingBases = baseCap - existingBases
	//			lastName, lastPos = findLastBinnedPos(answer[k])
	//
	//			for f := range genome {
	//				var currBases []dna.Base
	//				var currFa Fasta
	//
	//				if genome[f].Name == lastName {
	//					contigBases := len(genome[f].Seq) - lastPos
	//					if contigBases <= remainingBases { //TODO: this wouldn't ever happen?
	//						currBases = append(currBases, genome[f].Seq[lastPos+1:])
	//						currFa.Name = lastName
	//						currFa.Seq = currBases
	//						answer[k] = append(answer[k], currFa)
	//					} else if contigBases > remainingBases {
	//						currBases = append(currBases, genome[f].Seq[lastPos+1:lastPos+remainingBases])
	//						currFa.Name = lastName
	//						currFa.Seq = currBases
	//						answer[k] = append(answer[k], currFa)
	//					}
	//				}
	//				lenOfContig := len(genome[f].Seq)
	//				if lenOfContig > remainingBases {
	//
	//					currFa.Name = genome[f].Name + "0-" + fileio.IntToString(baseCap-1) //-1 since bases are zero-based
	//					currBases = append(currBases, genome[f].Seq[0:remainingBases-1])
	//					//TODO: make the above line into helper fucntion
	//					currFa.Seq = currBases
	//				} else if lenOfContig == remainingBases {
	//
	//				} else if lenOfContig < remainingBases {
	//
	//				}
	//			}
	//		}
	//	}
	//}

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
