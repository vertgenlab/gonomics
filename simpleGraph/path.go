package simpleGraph

import (
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/fastq"
	"log"
	"strings"
	"sync"
)

func PathToSeq(alignPath []uint32, samfile *sam.SamAln, gg *SimpleGraph) []dna.Base {
	var answer []dna.Base
	if alignPath == nil {
		answer = UnMappedRead(len(samfile.Seq))
	} else {
		for i := 0; i < len(alignPath); i++ {
			answer = append(answer, gg.Nodes[alignPath[i]].Seq...)
		}
	}
	return answer
}

func SamToPath(aln *sam.SamAln) []uint32 {

	words := strings.Split(aln.Extra, "\t")
	if len(words) < 2 {
		//log.Fatalf("Error: input sam might not contain path information...")
		//log.Printf("Possible unmapped read:%s\n", sam.SamAlnToString(aln))
		return nil
	} else {
		return StringToPath(words[1])
	}
	//log.Printf("%d\n", StringToPath(words[1]))

}

/*
func betterSamToPath(aln *sam.SamAln) []uint32 {
	words := strings.Split(aln.Extra, "\t")
	answer := make([]uint32, 0)
	if len(words) < 2 || sam.IsUnmapped(aln) {
		//log.Fatalf("Error: input sam might not contain path information...")
		//log.Printf("Possible unmapped read:%s\n", sam.SamAlnToString(aln))
		return nil
	}
	words[1] = words[1][5:]
	words = strings.Split(words[1], ":")
	for i := 0; i < len(words); i++ {
		answer = append(answer, common.StringToUint32(words[i]))
	}
	return answer
	//return StringToPath(words[1])
}*/

func StringToPath(allPaths string) []uint32 {
	words := strings.Split(allPaths[5:], ":")
	answer := make([]uint32, len(words))
	if strings.Compare(words[0], "-1") == 0 {
		return nil
	} else {
		for i := 0; i < len(words); i++ {
			answer[i] = common.StringToUint32(words[i])
		}
	}
	return answer
}

func SamChanView(incomingSams <-chan *sam.SamAln, gg *SimpleGraph, wg *sync.WaitGroup) {
	var yes, no, numReads int = 0, 0, 0
	log.SetFlags(log.Ldate | log.Ltime)
	for alignedRead := range incomingSams {
		numReads++
		//log.Printf("path=%v\n", SamToPath(alignedRead))

		log.Printf("%s\n", ViewGraphAlignment(alignedRead, gg))
		if CheckAlignment(alignedRead) {
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

func UnMappedRead(length int) []dna.Base {
	answer := make([]dna.Base, length)
	//var answer *fastq.Fastq = {Name: "UnMapped_Read", Seq: make([]dna.Base, length)
	for i := 0; i < len(answer); i++ {
		answer[i] = dna.N
	}
	return answer
}

func ViewGraphAlignment(samLine *sam.SamAln, genome *SimpleGraph) string {
	if SamToPath(samLine) == nil {
		return fmt.Sprintf("Unmapped Alignment:\n%s\n", sam.ModifySamToString(samLine, false, false, true, false, true, false, false, false, false, false, true))
	} else {

		var seqOne, seqTwo bytes.Buffer
		var operations []*cigar.Cigar = samLine.Cigar
		var i int64 = samLine.Pos - 1
		var j int64 = getStartRead(samLine)
		//var k int64
		var count int64
		//words := strings.Split(samLine.RName, "_")
		//log.Printf("path=%v\n", SamToPath(samLine))
		var alpha []dna.Base = PathToSeq(SamToPath(samLine), samLine, genome)
		var beta []dna.Base = samLine.Seq
		//var lineLength int = 50
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
		var lineLength int64 = 50
		var k, pos int64
		var prettySeq string = ""
		pos = addStartChrPos(samLine) + samLine.Pos
		for k = 0; k < int64(len(seqOne.String())); k += lineLength {

			if k+lineLength > int64(len(seqOne.String())) && k+lineLength > int64(len(seqTwo.String())) {

				prettySeq += fmt.Sprintf("%s:\t[%d-%d]\n", samLine.RName, k+pos, k+lineLength+pos) + fmt.Sprintf("%s\n", seqOne.String()[k:]) + fmt.Sprintf("%s\n", seqTwo.String()[k:])
			} else {
				prettySeq += fmt.Sprintf("%s:\t[%d-%d]\n", samLine.RName, k+pos, k+lineLength+pos) + fmt.Sprintf("%s\n", seqOne.String()[k:k+lineLength]) + fmt.Sprintf("%s\n", seqTwo.String()[k:k+lineLength])
			}
		}
		return fmt.Sprintf("%s\n%s", ModifySamToString(samLine, false, false, false, false, true, false, false, false, false, false, true), prettySeq)
	}
}

func addStartChrPos(samfile *sam.SamAln) int64 {
	var answer int64 = 0
	if strings.Contains(samfile.Extra, "XO:i:") {
		words := strings.Split(samfile.Extra, "\t")
		answer = common.StringToInt64(words[2][5:])
	}
	return answer
}

func ModifySamToString(aln *sam.SamAln, samflag bool, rname bool, pos bool, mapq bool, cig bool, rnext bool, pnext bool, tlen bool, seq bool, qual bool, extra bool) string {
	var answer string = fmt.Sprintf("%s\n", aln.QName)
	if samflag {
		answer += fmt.Sprintf("%d\n", aln.Flag)
	}
	if rname {
		answer += fmt.Sprintf("%s\t", aln.RName)
	}
	if pos {
		aln.Pos += addStartChrPos(aln)
		answer += fmt.Sprintf("%d\t", aln.Pos)
	}
	if mapq {
		answer += fmt.Sprintf("%d\t", aln.MapQ)
	}
	if cig {
		answer += fmt.Sprintf("%s\t", cigar.ToString(aln.Cigar))
	}
	if rnext {
		answer += fmt.Sprintf("%s\t", aln.RNext)
	}
	if pnext {
		answer += fmt.Sprintf("%d\t", aln.PNext)
	}
	if tlen {
		answer += fmt.Sprintf("%d\t", aln.TLen)
	}
	if seq {
		answer += fmt.Sprintf("%s\t", dna.BasesToString(aln.Seq))
	}
	if qual {
		answer += fmt.Sprintf("%s\t", string(aln.Qual))
	}
	if extra {
		words := strings.Split(aln.Extra, "\t")
		for _, text := range words[:len(words)-1] {
			if strings.Contains(text, "GP:Z:") {

				answer += fmt.Sprintf("GP:Z:\n%s", pathPrettyString(text[5:]))
			} else {
				answer += fmt.Sprintf("%s\t", text)
			}
		}
	}
	return answer
}

func pathPrettyString(graphPath string) string {

	var s string = ""

	words := strings.Split(graphPath, ":")
	//log.Printf("%v\n", words)
	var i int
	var j int
	for i = 0; i < len(words); i += 8 {
		var line string = ""
		if i+8 > len(words) {
			line += fmt.Sprintf("%s", words[i])
			for j = i + 1; j < len(words)-1; j++ {
				line += fmt.Sprintf(":%s", words[j])
			}
			s += fmt.Sprintf("%s\n", line)
		} else {
			line += fmt.Sprintf("%s", words[i])
			for j = i + 1; j < i+8; j++ {
				line += fmt.Sprintf(":%s", words[j])
			}
			s += fmt.Sprintf("%s\n", line)
		}
		//s += fmt.Sprintf("\n")
	}

	return s
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
