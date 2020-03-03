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
)

func PathToSeq(alignPath []uint32, gg *SimpleGraph, start int64) []dna.Base {
	answer := make([]dna.Base, 0)
	answer = append(answer, gg.Nodes[alignPath[0]].Seq[start:]...)
	if len(alignPath) > 1 {
		for i := 1; i < len(alignPath); i++ {
			answer = append(answer, gg.Nodes[alignPath[i]].Seq...)
		}
	}
	return answer
}

func ViewGraphAignment(samLine *sam.SamAln, genome *SimpleGraph) string {
	//var maxI int64
	//var operations []*Cigar
	var seqOne, seqTwo bytes.Buffer

	var operations []*cigar.Cigar = samLine.Cigar
	var i int64 = 0
	var j int64 = 0
	var count int64
	//words := strings.Split(samLine.RName, "_")
	var alpha []dna.Base = PathToSeq(SamToPath(samLine), genome, samLine.Pos-1)
	var beta []dna.Base = samLine.Seq

	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case 'M':
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
		log.Fatalf("Error: input sam might not contain path information...")
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
