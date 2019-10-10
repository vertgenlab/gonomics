package simpleGraph

import (
    //"fmt"
    "math/rand"
    "github.com/vertgenlab/gonomics/fasta"
    //"github.com/vertgenlab/gonomics/dna"
    "strconv"
)
func randIntInRange(x int, y int) int {
    return int(rand.Float64() * float64(y - x)) + x
}

func RandomSeqGenerator(genome []*fasta.Fasta, readLength int, readNumber int) []*fasta.Fasta {
	var answer []*fasta.Fasta
	var curr *fasta.Fasta
	var start int
	//var startPos string
	//var endPos string
	var randStart int 
	var readName string
	for i := 0; i < readNumber; i++ {
		randStart = randIntInRange(0, len(genome))
		curr = genome[randStart]
		start = randIntInRange(0, len(genome[randStart].Seq)-readLength)
		readName = curr.Name + "_" + strconv.Itoa(start) + "_" + strconv.Itoa(start+readLength)
		
		//for j := 0; j < len(curr.Seq); j++ {
		//	if curr.Seq[j] == dna.N {
		//		i++
		//	}
		//}
		answer = append(answer, &fasta.Fasta{Name: readName, Seq: curr.Seq[start:start+readLength]})
	}
	answer = fasta.RemoveGaps(answer)
	fasta.Write("simReads.fasta", answer)
	return answer
}
