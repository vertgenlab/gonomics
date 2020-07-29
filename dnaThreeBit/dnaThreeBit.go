package dnaThreeBit

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

type ThreeBit struct {
	Seq []uint64
	Len int
}

// This could have been uint8, but with that I found lots
// of casting between uint64 and uint8.  I don't think many
// ThreeBitBases will be sitting around for a long time by themselves
// so I don't think the extra memory will be noticed.
type ThreeBitBase uint64

const (
	A          ThreeBitBase = 0
	C          ThreeBitBase = 1
	G          ThreeBitBase = 2
	T          ThreeBitBase = 3
	N          ThreeBitBase = 4
	PaddingOne ThreeBitBase = 5
	PaddingTwo ThreeBitBase = 6
)

//basesTo3bit
//3bitToBases
//stringTo3Bit
//3BitToString
//ReadFasta
//WriteFasta
//Read3bit
//Write3bit

func ReadFromFasta(filename string) []*ThreeBit {
        var line string
        var currSeq []dna.Base
        var answer []*Fasta
        var seqIdx int64 = -1
        var doneReading bool = false

        file := fileio.EasyOpen(filename)
        defer file.Close()

        for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
                if strings.HasPrefix(line, ">") {
                        tmp := Fasta{Name: line[1:len(line)]}
                        answer = append(answer, &tmp)
                        seqIdx++
                } else {
                        currSeq = dna.StringToBases(line)
                        answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq...)
                }
        }
        return answer
}

func NextSeqFasta(reader *fileio.EasyReader) *ThreeBit {
        var line, seqSoFar string
        var err error
        var nextBytes []byte
        var answer []dna.Base
        for nextBytes, err = reader.Peek(1); err == nil && nextBytes[0] != '>'; nextBytes, err = reader.Peek(1) {
                line, _ = fileio.EasyNextLine(reader)
                answer = append(answer, dna.StringToBases(line)...)
        }
        return answer
}

func BasesToUint64(seq []dna.Base, start int, end int, padding ThreeBitBase) uint64 {
	if end-start > 21 || start >= end {
		log.Fatalf("Error: when converting to ThreeBit. start=%d end=%d\n", start, end)
	}
	var idx int
	var answer uint64 = uint64(seq[start]) << 1
	for idx = start + 1; idx < end; idx++ {
		answer = answer << 3
		answer = answer | (uint64(seq[idx]) << 1)
	}
	for ; idx < start+21; idx++ {
		answer = answer << 3
		answer = answer | (uint64(padding) << 1)
	}
	return answer
}

// basesToUint64WithOffset places padding at the beginning of the sequence (as well as the end, which is normal)
// so that sequences may be "in register" with each other so that an xor can quickly compare them for equality
func basesToUint64WithOffset(seq []dna.Base, start int, end int, padding ThreeBitBase, offset int) uint64 {
        if end-start+offset > 21 || start >= end {
                log.Fatalf("Error: when converting to ThreeBit. start=%d end=%d\n", start, end)
        }
        var idx int
	var answer uint64 = 0
	for idx = 0; idx < offset; idx++ {
		answer = answer << 3 // not needed the first time through the loop, but does not hurt
		answer = answer | (uint64(padding) << 1)
	}
        for idx = start; idx < end; idx++ {
                answer = answer << 3
                answer = answer | (uint64(seq[idx]) << 1)
        }
        for ; idx < start+21; idx++ {
                answer = answer << 3
                answer = answer | (uint64(padding) << 1)
        }
        return answer
}

func GetThreeBitBase(fragment *ThreeBit, pos int) ThreeBitBase {
        if pos < 0 || pos >= fragment.Len {
                log.Fatalf("Error: asked for base at position:%d for a sequence with length:%d\n", pos, fragment.Len)
        }
        var upos = uint(pos)
        const lastBase uint64 = 7 // right-most three bits are one
        var idx uint = upos / 21
        var remainder uint = upos % 21
        var shift uint = 64 - 3*(remainder+1)
        return (fragment.Seq[idx] >> shift) & lastBase
}

func GetBase(fragment *ThreeBit, pos int) dna.Base {
	return dna.Base(GetThreeBitBase(fragment, pos))
}

func NewThreeBit(inSeq []dna.Base, padding ThreeBitBase) *ThreeBit {
	var sliceLenNeeded int = (len(inSeq) + 20) / 21
	var start, end int = 0, 0
	answer := ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)}
	for i := 0; i < sliceLenNeeded; i++ {
		start = i * 21
		end = common.Min(start+21, len(inSeq))
		answer.Seq[i] = BasesToUint64(inSeq, start, end, padding)
	}
	return &answer
}

func newThreeBitWithOffset(inSeq []dna.Base, padding ThreeBitBase, offset int) *ThreeBit {
        var sliceLenNeeded int = (len(inSeq) + offset + 20) / 21
        var start, end int = 0, 0
        answer := ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)+offset}
        for i := 0; i < sliceLenNeeded; i++ {
                start = i * 21
                end = common.Min(start+21, len(inSeq))
                answer.Seq[i] = BasesToUint64(inSeq, start, end, padding)
        }
        return &answer
}
