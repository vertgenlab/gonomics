package fasta

import (
	"math/rand"
	"strconv"

	"github.com/vertgenlab/gonomics/dna"
)

func randIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandomSeqGenerator(genome []Fasta, readLength int, readNumber int) []Fasta {
	var answer []Fasta
	//var curr *Fasta
	var start int
	//var startPos string
	//var endPos string
	var chrom int
	var readName string
	for i := 0; i < readNumber; {
		chrom = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chrom].Seq)-readLength)
		readName = genome[chrom].Name + "_" + strconv.Itoa(start) + "_" + strconv.Itoa(start+readLength)

		if dna.CountBase(genome[chrom].Seq, dna.N) == 0 {
			answer = append(answer, Fasta{Name: readName, Seq: genome[chrom].Seq[start : start+readLength]})
			i++
		} else {

		}
	}
	//answer = fasta.RemoveGaps(answer)
	//fasta.Write("simReads.fasta", answer)
	return answer
}
