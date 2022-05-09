package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"io"
)

func (s Sam) GetChrom() string {
	return s.RName
}

func (s Sam) GetChromStart() int {
	return int(s.Pos - 1)
}

func (s Sam) GetChromEnd() int {
	var runLength int = 0
	if s.Cigar[0].Op == '*' {
		return s.GetChromStart()
	}
	for i := 0; i < len(s.Cigar); i++ {
		if cigar.ConsumesReference(s.Cigar[i].Op) {
			runLength += int(s.Cigar[i].RunLength)
		}
	}
	return s.GetChromStart() + runLength
}

func (s Sam) UpdateLift(c string, start int, end int) {
	s.RName = c
	s.Pos = uint32(start) + 1
}

func (s Sam) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, s)
}
