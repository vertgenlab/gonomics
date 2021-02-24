package simpleGraph

import (
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

func PathToSeq(p giraf.Path, genome *SimpleGraph) []dna.Base {
	if len(p.Nodes) == 1 {
		answer := make([]dna.Base, p.TEnd-p.TStart)
		copy(answer, genome.Nodes[p.Nodes[0]].Seq[p.TStart:p.TEnd])
		return answer
	} else {
		answer := make([]dna.Base, len(genome.Nodes[p.Nodes[0]].Seq)-p.TStart)
		copy(answer, genome.Nodes[p.Nodes[0]].Seq[p.TStart:])
		for i := 1; i < len(p.Nodes)-1; i++ {
			answer = append(answer, genome.Nodes[p.Nodes[i]].Seq...)
		}
		answer = append(answer, genome.Nodes[p.Nodes[len(p.Nodes)-1]].Seq[:p.TEnd]...)
		return answer
	}
}

func ViewGraphAlignment(g *giraf.Giraf, genome *SimpleGraph) string {

	var seqOne, seqTwo bytes.Buffer
	var i int = g.Path.TStart
	var j int = g.QStart

	var count int
	var alpha []dna.Base = PathToSeq(g.Path, genome)
	var beta []dna.Base = g.Seq
	for _, operation := range g.Cigar {
		for count = 0; count < int(operation.RunLen); count++ {
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
	} /*
		var lineLength int64 = 50
		var k, pos int
		var prettySeq string = ""
		pos = addStartChrPos(samLine) + samLine.Pos
		for k = 0; k < int64(len(seqOne.String())); k += lineLength {

			if k+lineLength > int64(len(seqOne.String())) && k+lineLength > int64(len(seqTwo.String())) {

				prettySeq += fmt.Sprintf("%s:\t[%d-%d]\n", samLine.RName, k+pos, k+lineLength+pos) + fmt.Sprintf("%s\n", seqOne.String()[k:]) + fmt.Sprintf("%s\n", seqTwo.String()[k:])
			} else {
				prettySeq += fmt.Sprintf("%s:\t[%d-%d]\n", samLine.RName, k+pos, k+lineLength+pos) + fmt.Sprintf("%s\n", seqOne.String()[k:k+lineLength]) + fmt.Sprintf("%s\n", seqTwo.String()[k:k+lineLength])
			}
		}*/
	//return fmt.Sprintf("%s\n%s", ModifySamToString(samLine, false, true, true, false, true, false, false, false, false, false, true), prettySeq)
	return fmt.Sprintf("%s\n%s", seqOne.String(), seqTwo.String())
}

/*
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
}*/

func addStartChrPos(samfile *sam.SamAln) int64 {
	var answer int64 = 0
	if strings.Contains(samfile.Extra, "XO:i:") {
		words := strings.Split(samfile.Extra, "\t")
		answer = common.StringToInt64(words[2][5:])
	}
	return answer
}

func ModifySamToString(aln *sam.SamAln, samflag bool, rname bool, pos bool, mapq bool, cig bool, rnext bool, pnext bool, tlen bool, seq bool, qual bool, extra bool) string {
	var answer string = fmt.Sprintf("%s\n\n", aln.QName)
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
	var i, j int
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
	}
	return s
}

func AddPath(allPaths []uint32, newPath uint32) []uint32 {
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
		currPaths = AddPath(currPaths, newPaths[0])
		currPaths = append(currPaths, newPaths[1:]...)
		return currPaths
	}
}

func reversePath(alpha []uint32) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func PathToString(allPaths []uint32) string {
	var s string = ""
	if allPaths == nil {
		return s
	} else {
		s += fmt.Sprint(allPaths[0])
		if len(allPaths) > 1 {
			for i := 1; i < len(allPaths); i++ {
				s += ":" + fmt.Sprint(allPaths[i])
			}
		}
	}
	return s
}

func getSeedPath(seed *SeedDev) []uint32 {
	var path []uint32 = []uint32{seed.TargetId}
	if seed.NextPart == nil {
		return path
	} else {
		path = append(path, getSeedPath(seed.NextPart)...)
	}
	return path
}

func getStartRead(aln *sam.SamAln) int64 {
	var alignedPos int = 0
	if aln.Cigar[0].Op == 'S' {
		alignedPos += aln.Cigar[0].RunLength
	}
	return int64(alignedPos)
}
