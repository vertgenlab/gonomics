package genomeGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
)

// TODO: what about neg strand?
func perfectMatchBig(read fastq.FastqBig, scoreMatrix [][]int64) int64 {
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

func scoreSeedFastqBig(seed *Seed, read fastq.FastqBig, scoreMatrix [][]int64) int64 {
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

func scoreSeedPart(seed *Seed, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreSeed(seed *Seed, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read, scoreMatrix)
	}
	return score
}

/*func NodesHeader(ref []*Node) *sam.Header {
	var header sam.Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s_%d\tLN:%d", ref[i].Name, ref[i].Id, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: len(ref[i].Seq)})
	}
	return &header
}*/
