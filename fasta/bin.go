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

func fewerContigsThanBins(fa []Fasta, binNum int) map[int][]Fasta {
	var whichBin int
	var remainingBases int
	var answer map[int][]Fasta
	totalBases := 0
	for c := range fa {
		totalBases = totalBases + len(fa[c].Seq)
	}
	//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with baseCap,
	//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
	//	could do a check to see if remainder is super large and then handle as necessary
	equalBinNum := binNum - 1
	baseCap := totalBases / equalBinNum
	leftOvers := totalBases % equalBinNum

	for k := range answer {
		existingBases := 0
		if answer[k] != nil {
			existingBases = calcNumBasesInBin(answer[k])
			if existingBases < baseCap {
				whichBin = k
				remainingBases = baseCap - existingBases
			}
		}
	}
	for f := range fa {
		var currFa Fasta
		var currBases []dna.Base
		lenOfContig := len(fa[f].Seq)
		if lenOfContig > baseCap {

			currFa.Name = fa[f].Name + "0-" + fileio.IntToString(baseCap-1) //-1 since bases are zero-based
			currBases = append(currBases, fa[f].Seq[0:baseCap-1])
			//TODO: make the above line into helper fucntion
			currFa.Seq = currBases
			for k := range answer {
				existingBases := 0
				if answer[k] != nil {
					existingBases = calcNumBasesInBin(answer[k])
				}
				if existingBases < baseCap-1 {
					answer[k] = append(answer[k], currFa)
				}
			}
		} else if lenOfContig == baseCap {

		} else if lenOfContig < baseCap {

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
