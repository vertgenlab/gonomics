package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"math/rand"
	"strconv"
)

func randIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandomReads(genome []*Node, readLength int, numReads int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var readName string
	var strand bool
	var qual []rune
	var seq []dna.Base
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, int64(start), int64(start+readLength)) == 0 {
			readName = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			qual = generateFakeQual(readLength)
			seq = genome[chromIdx].Seq[start : start+readLength]
			if !strand {
				dna.ReverseComplement(seq)
			}
			answer[i] = &fastq.Fastq{Name: readName, Seq: seq, Qual: qual}
			i++
		}
	}
	return answer
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

		if dna.CountBase(genome[chrom].Seq[start:start+readLength], dna.N) == 0 {
			qual = generateQual(genome[chrom].Seq[start : start+readLength])
			seq = genome[chrom].Seq[start : start+readLength]
			dna.AllToUpper(seq)
			answer = append(answer, &fastq.Fastq{Name: readName, Seq: seq, Qual: qual})
			i++
		}
	}
	return answer
}

func generateFakeQual(length int) []rune {
	var answer []rune = make([]rune, length)
	for i := 0; i < length; i++ {
		answer[i] = 'J'
	}
	return answer
}

func generateQual(bases []dna.Base) []rune {
	var ans []rune
	for i := 0; i < len(bases); i++ {
		ans = append(ans, 'J')
	}
	return ans
}
