package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

//TODO: what about neg strand?
func perfectMatchBig(read *fastq.FastqBig, scoreMatrix [][]int64) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func scoreSeedSeq(seq []dna.Base, start uint32, end uint32, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := start; i < end; i++ {
		score += scoreMatrix[seq[i]][seq[i]]
	}
	return score
}

func scoreSeedFastqBig(seed *SeedDev, read *fastq.FastqBig, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		if seed.PosStrand {
			score += scoreMatrix[read.Seq[i]][read.Seq[i]]
		} else {
			score += scoreMatrix[read.SeqRc[i]][read.SeqRc[i]]
		}
	}
	return score
}

func scoreSeedPart(seed *SeedDev, read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreSeed(seed *SeedDev, read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read, scoreMatrix)
	}
	return score
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

func SoftClipBases(front int, lengthOfRead int, cig []cigar.ByteCigar) []cigar.ByteCigar {
	var runLen int = cigar.QueryRunLen(cig)
	if runLen < lengthOfRead {
		answer := make([]cigar.ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, cigar.ByteCigar{RunLen: uint32(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+cigar.QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, cigar.ByteCigar{RunLen: uint32(lengthOfRead - front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}

}

func NewSwMatrix(size int) *MatrixAln {
	sw := MatrixAln{}
	sw.m, sw.trace = MatrixSetup(size)
	return &sw
}

func MatrixSetup(size int) ([][]int64, [][]byte) {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return m, trace
}

//perfect match
func perfectMatch(read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func NodesHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s_%d\tLN:%d", ref[i].Name, ref[i].Id, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})
	}
	return &header
}

func ChromAndPosToNumber(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func dnaToNumber(seq []dna.Base, start int, end int) uint64 {
	var answer uint64 = uint64(seq[start])
	for i := start + 1; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(seq[i])
	}
	return answer
}

func numberToChromAndPos(code uint64) (int64, int64) {
	var rightSideOnes uint64 = 4294967295
	var leftSideOnes uint64 = rightSideOnes << 32
	var chromIdx uint64 = code & leftSideOnes
	chromIdx = chromIdx >> 32
	var pos uint64 = code & rightSideOnes
	return int64(chromIdx), int64(pos)
}

func ViewMatrix(m [][]int64) string {
	var message string = ""
	message += fmt.Sprintf("\t\t %d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t %d\n", m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3])
	return message
}

func CheckAlignment(aln *sam.SamAln, genome *SimpleGraph) bool {
	var answer bool = false
	samPath := SamToPath(aln)
	if samPath == nil {
		return false
	}
	qName := strings.Split(aln.QName, "_")
	if strings.Compare(genome.Nodes[common.StringToUint32(qName[0])].Name, aln.RName) == 0 && aln.Pos == common.StringToInt64(qName[1]) {
		return true
	}
	return answer
}

func CheckAnswers(query []*sam.SamAln, genome *SimpleGraph) {
	var yes, no int64 = 0, 0
	for i := 0; i < len(query); i++ {
		if CheckAlignment(query[i], genome) {
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
