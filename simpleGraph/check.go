package simpleGraph

import (
	//"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	//"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/dna"
	"log"
	//"os"
	"strings"
	//"sync"
	//"time"
)

func checkAlignment(aln *sam.SamAln) bool {
	var answer bool = false
	words := strings.Split(aln.QName, "_")
	var ref int64
	//var query int64
	var blastScore int64
	alignedPos := common.StringToInt64(words[1])
	if alignedPos <= ref-10 || (alignedPos > ref+100) {
		words = strings.Split(aln.Extra, "\t")
		blastScore = common.StringToInt64(words[0][5:])
		if blastScore > 5000 {
			answer = true
		}
		log.Printf("\t%s\t%s\t%d\t%s\t%s\n", cigar.ToString(aln.Cigar), aln.QName, aln.Pos, aln.RName, aln.Extra)
		//cigar.UpdateIndexSlice(ref, 0, aln.Cigar)
		//if simSam.Pos+ref == common.StringToInt64(words[2]) {
		//	answer = true
		//}
	}
	return answer
}

func CheckAnswers(query []*sam.SamAln) {
	var yes, no int64 = 0, 0
	for i := 0; i < len(query); i++ {
		if checkAlignment(query[i]) {
			yes++
			//log.Printf(sam.SamAlnToString(query[i]))
		} else {
			no++
			//log.Printf("This did not map:\n%s\n", sam.SamAlnToString(query[i]))
		}
	}
	log.Printf("Total number of reads aligned: %d...", len(query))
	log.Printf("Number of reads correctly aligned: %d...\n", yes)
	log.Printf("Number of reads mismapped: %d...\n", no)
}

/*
func WriteReadsToFile(genome []*Node, readLength int, numReads int) {

	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 1

	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("gasAcu1.fa")

	log.Printf("Simulating %d reads...\n", numberOfReads)
	simReads := RandomReads(genome, readLength, numberOfReads, mutations, true)
	simReads := GenomeDiversity(genome.Nodes, readLength, numberOfReads)
	fastq.Write("genomeDiversity.fq", simReads)
	log.Printf("Finished writing to file %d reads...\n", numberOfReads)
}*/
