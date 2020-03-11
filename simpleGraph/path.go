package simpleGraph

import (
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
	"sync"
)

func PathToSeq(alignPath []uint32, gg *SimpleGraph) []dna.Base {
	var answer []dna.Base
	for i := 0; i < len(alignPath); i++ {
		answer = append(answer, gg.Nodes[alignPath[i]].Seq...)
	}
	return answer
}

func SamChanView(incomingSams <-chan *sam.SamAln, gg *SimpleGraph, wg *sync.WaitGroup) {
	var yes, no, numReads int = 0, 0, 0
	for alignedRead := range incomingSams {
		numReads++
		//log.Printf("path=%v\n", SamToPath(alignedRead))

		log.Printf("%s\n", ViewGraphAignment(alignedRead, gg))
		if checkAlignment(alignedRead) {
			yes++
		} else {
			no++
		}
	}
	log.Printf("Total number of reads aligned: %d...", numReads)
	log.Printf("Number of reads correctly aligned: %d...\n", yes)
	log.Printf("Number of reads mismapped: %d...\n", no)
	wg.Done()
}

func ViewGraphAignment(samLine *sam.SamAln, genome *SimpleGraph) string {
	var seqOne, seqTwo bytes.Buffer

	var operations []*cigar.Cigar = samLine.Cigar
	var i int64 = samLine.Pos - 1
	var j int64 = getStartRead(samLine)
	var count int64
	//words := strings.Split(samLine.RName, "_")
	//log.Printf("path=%v\n", SamToPath(samLine))
	var alpha []dna.Base = PathToSeq(SamToPath(samLine), genome)
	var beta []dna.Base = samLine.Seq

	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case 'M':
				if i >= int64(len(alpha)) {
					log.Printf("%s\n", sam.SamAlnToString(samLine))
					log.Printf("target:%s\n", dna.BasesToString(alpha))
					log.Printf("query:%s\n", dna.BasesToString(beta))
				}
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
			case 'I':
				seqOne.WriteRune('-')
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				j++
			case 'D':
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune('-')
				i++
			case 'S':
				seqOne.WriteRune('-')
				seqTwo.WriteRune('-')
			}
		}
	}
	return fmt.Sprintf("%s:%s:%s\n%s\n%s\n", samLine.QName, cigar.ToString(operations), samLine.Extra, seqOne.String(), seqTwo.String())
}

func SamToPath(aln *sam.SamAln) []uint32 {

	words := strings.Split(aln.Extra, "\t")
	if len(words) < 2 {
		//log.Fatalf("Error: input sam might not contain path information...")
		log.Printf("Possible unmapped read:%s\n", sam.SamAlnToString(aln))
		return nil
	}
	//log.Printf("%d\n", StringToPath(words[1]))
	return StringToPath(words[1])
}

func AddPath(newPath uint32, allPaths []uint32) []uint32 {
	if len(allPaths) == 0 {
		allPaths = append(allPaths, newPath)
	} else if allPaths[len(allPaths)-1] == newPath {
		return allPaths
	} else {
		allPaths = append(allPaths, newPath)
	}
	return allPaths
}

func CatPaths(currPaths []uint32, newPaths []uint32) []uint32 {
	if len(newPaths) == 0 {
		return currPaths
	} else if len(currPaths) == 0 {
		return newPaths
	} else {
		currPaths = AddPath(newPaths[0], currPaths)
		currPaths = append(currPaths, newPaths[1:]...)
		return currPaths
	}
}

func reversePath(alpha []uint32) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func PathToString(allPaths []uint32, gg *SimpleGraph) string {
	var s string = ""
	//fmt.Printf("length of paths %d\n", len(allPaths))
	if allPaths == nil {
		return s
	} else {
		s += fmt.Sprint(gg.Nodes[allPaths[0]].Id)
		if len(allPaths) > 1 {
			for i := 1; i < len(allPaths); i++ {
				s += ":" + fmt.Sprint(gg.Nodes[allPaths[i]].Id)
			}
		}
	}
	return s
}

func StringToPath(allPaths string) []uint32 {
	words := strings.Split(allPaths, ":")
	answer := make([]uint32, len(words)-2)

	for i := 0; i < len(words)-2; i++ {
		answer[i] = common.StringToUint32(words[i+2])
	}
	return answer
}

//PathToString(CatPaths(CatPaths(reversePath(leftPath), getSeedPath(seeds[i])), rightPath), gg)
func getSeedPath(seed *SeedDev) []uint32 {
	var path []uint32 = []uint32{seed.TargetId}
	if seed.Next == nil {
		return path
	} else {
		path = append(path, getSeedPath(seed.Next)...)
	}
	return path
}

func getStartRead(aln *sam.SamAln) int64 {
	var alignedPos int64 = 0
	if aln.Cigar[0].Op == 'S' {
		alignedPos += aln.Cigar[0].RunLength
	}
	return alignedPos
}
