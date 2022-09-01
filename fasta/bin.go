package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func BinFasta(genome []Fasta, binNum int) map[int][]Fasta {
	var totalBases, lastPos int
	var lastName string
	answer := make(map[int][]Fasta, binNum)
	for c := range genome {
		totalBases = totalBases + len(genome[c].Seq)
	}

	//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with baseCap,
	//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
	//	could do a check to see if remainder is super large and then handle as necessary
	equalBinNum := binNum - 1
	remainder := totalBases % binNum
	splitBases := totalBases - remainder //number of bases to split between equalBins
	baseCap := splitBases / equalBinNum

	for j := 0; j < binNum; j++ {
		for i := 0; i < len(genome); i++ {

			log.Print(answer)
			var remainingBases int
			var done int
			var currFa Fasta
			var currBases []dna.Base
			if j == binNum { //we have finished making all bins but the last bin which will hold the remainder of bases
				answer[j], done = fillLastBin(genome, i)
				i = done
			}

			value, exists := answer[j]

			log.Printf("evaluating answer[%d]", j)
			log.Print(value, exists)

			if j == 0 && i == 0 {
				currBases = appendRangeOfBases(genome[i].Seq, 0, baseCap)
				currFa.Name = genome[i].Name
				currFa.Seq = currBases
				answer[j] = append(answer[j], currFa)
				lastPos = baseCap
				continue
			} else if value != nil { //old bin

				log.Print("old bin new contig")

				filledBases := calcNumBasesInBin(answer[j])
				if filledBases < baseCap { //bin not full
					remainingBases = baseCap - filledBases

					log.Print(remainingBases)

					currBases = appendRangeOfBases(genome[i].Seq, 0, remainingBases)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
					lastPos = len(currBases)

					log.Print("bin not full", lastPos)

				} else {
					continue
				}
			} else if !exists || value == nil { //new bin

				log.Print("making new bin")
				log.Print(i)

				lastName = findLastBinnedName(answer[j-1])

				log.Print(lastPos, len(genome[i].Seq))

				if lastName == genome[i].Name && lastPos < len(genome[i].Seq) { //new bin and old contig

					log.Print("new bin old contig")

					remainingBases = len(genome[i].Seq) - lastPos
					currBases = appendRangeOfBases(genome[i].Seq, lastPos, lastPos+remainingBases)
					currFa.Name = genome[i].Name
					currFa.Seq = currBases
					answer[j] = append(answer[j], currFa)
					lastPos = lastPos + len(currBases)
				} else if lastName == genome[i].Name && lastPos == len(genome[i].Seq) { //new bin and new contig

					log.Print("new bin new contig")
					lastPos = 0

					//currBases = appendRangeOfBases(genome[i].Seq, 0, baseCap)
					//currFa.Name = genome[i].Name
					//currFa.Seq = currBases
					//answer[j] = append(answer[j], currFa)
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

	log.Print(answer)

	return answer
}

func calcNumBasesInBin(f []Fasta) int {
	totalBases := 0
	for i := range f {
		totalBases = totalBases + len(f[i].Seq)
	}
	return totalBases
}

func findLastBinnedName(f []Fasta) string {
	var lastFasta string

	for i := range f {
		lastFasta = f[i].Name
	}

	return lastFasta
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
			if b >= start && b <= stop-1 {
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
