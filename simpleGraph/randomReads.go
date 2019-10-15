package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"math/rand"
	"strconv"
)

func randIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandomFastqGen(genome []*fasta.Fasta, readLength int, readNumber int) []*fastq.Fastq {
	var answer []*fastq.Fastq
	//var curr *fastq.Fastq
	var start int
	//var startPos string
	//var endPos string
	var chrom int
	var readName string
	var qual []rune
	var seq []dna.Base
	for i := 0; i < readNumber; {
		chrom = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chrom].Seq)-readLength)
		readName = genome[chrom].Name + "_" + strconv.Itoa(start) + "_" + strconv.Itoa(start+readLength)

		if dna.CountBase(genome[chrom].Seq[start : start+readLength], dna.N) == 0 {
			qual = generateQual(genome[chrom].Seq[start : start+readLength])
			seq = genome[chrom].Seq[start : start+readLength]
			dna.AllToUpper(seq)
			answer = append(answer, &fastq.Fastq{Name: readName, Seq: seq, Qual: qual})
			i++
		}
	}
	return answer
}

func generateQual(bases []dna.Base) []rune {
	var ans []rune
	for i := 0 ; i < len(bases);i++ {
		ans = append(ans, 'J')
	}
	return ans
}

/*
func RandomSeqGenerator(genome []*fasta.Fasta, readLength int, readNumber int) []*fasta.Fasta {
	var answer []*fasta.Fasta
	//var curr *fasta.Fasta
	var start int
	//var startPos string
	//var endPos string
	var chrom int
	var seq []dna.Base
	var readName string
	for i := 0; i < readNumber; {
		chrom = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chrom].Seq)-readLength)
		readName = genome[chrom].Name + "_" + strconv.Itoa(start) + "_" + strconv.Itoa(start+readLength)

		if dna.CountBase(genome[chrom].Seq[start : start+readLength], dna.N) == 0 {
			seq = genome[chrom].Seq[start : start+readLength]
			dna.AllToUpper(seq)
			answer = append(answer, &fasta.Fasta{Name: readName, Seq: seq})
			i++
		}
	}
	//answer = fasta.RemoveGaps(answer)
	//fasta.Write("simReads.fasta", answer)
	return answer
}*/


