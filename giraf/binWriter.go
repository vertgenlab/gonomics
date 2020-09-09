package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
)

func CompressGiraf(g Giraf) BinGiraf {
	var answer BinGiraf
	answer.qNameLen = uint8(len(g.QName))
	answer.qName = g.QName
	answer.flag = g.Flag
	answer.tStart = uint32(g.Path.TStart)
	answer.tEnd = uint32(g.Path.TEnd)
	answer.pathLen = uint16(len(g.Path.Nodes))
	answer.path = g.Path.Nodes
	answer.numCigarOps = uint16(len(g.Cigar))
	answer.byteCigar = g.Cigar
	answer.fancySeq = getFancySeq(g.Seq, g.Cigar)
	answer.alnScore = int64(g.AlnScore)
	answer.mapQ = g.MapQ
	answer.qual = encodeQual(g.Qual)
	answer.numQualOps = uint16(len(answer.qual))
	answer.notes = encodeNotes(g.Notes)
	//TODO hasNamePrefix
	//TODO blocksize
	return answer
}

func getFancySeq(seq []dna.Base, cigar []cigar.ByteCigar) dnaThreeBit.ThreeBit {
	var answer []dna.Base
	var seqIdx int
	for _, val := range cigar {
		if val.Op == 'S' || val.Op == 'X' || val.Op == 'I' {
			answer = append(answer, seq[seqIdx:seqIdx+int(val.RunLen)]...)
		}
		seqIdx += int(val.RunLen)
	}
	return *dnaThreeBit.NewThreeBit(answer, dnaThreeBit.A)
}

func encodeQual(q []uint8) []cigar.ByteCigar {
	answer := make([]cigar.ByteCigar, 0, len(q))
	var curr cigar.ByteCigar
	curr.Op = q[0]
	for i := 0; i < len(q); i++ {
		if q[i] != curr.Op && curr.RunLen != 0 {
			answer = append(answer, curr)
			curr.RunLen = 0
			curr.Op = q[i]
		}
		curr.RunLen++
	}

	if curr.RunLen != 0 {
		answer = append(answer, curr)
	}

	return answer
}

func encodeNotes(n []Note) []BinNote {
	answer := make([]BinNote, len(n))
	for idx, note := range n {
		answer[idx] = BinNote{
			tagType: note.Type,
			data:    []byte(note.Value)}
		copy(answer[idx].tag[:], note.Tag[:2])
	}
	return answer
}
