package sam

import (
	"io"

	"github.com/vertgenlab/gonomics/cigar"
)

func (s Sam) GetChrom() string {
	return s.RName
}

func (s Sam) GetChromStart() int {
	return int(s.Pos - 1)
}

func (s Sam) GetChromEnd() int {
	var runLength int = 0
	if s.Cigar == nil || s.Cigar[0].Op == '*' {
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

// Calculate the estimated initial capacity
func estimateInitialCapacity(s Sam) int {
	estimatedLength := 0

	// Fixed text parts
	fixedParts := []string{
		"QName: ", ", Flag: ", ", MapQ: ", ", RName: ", ", Pos: ", ", Cigar: ",
		", RNext: ", ", PNext: ", ", TLen: ", ", Seq: ", ", Qual: ", ", Extra: ",
	}
	for _, part := range fixedParts {
		estimatedLength += len(part)
	}

	// Dynamic parts
	estimatedLength += len(s.QName)
	estimatedLength += 5 // for Flag, assuming max value 65535
	estimatedLength += 3 // for MapQ, assuming max value 255
	estimatedLength += len(s.RName)
	estimatedLength += 10                // for Pos, assuming max value 4294967295
	estimatedLength += len(s.Cigar) * 10 // approximate size for each Cigar element
	estimatedLength += len(s.RNext)
	estimatedLength += 10             // for PNext, assuming max value 4294967295
	estimatedLength += 11             // for TLen, assuming range -2147483648 to 2147483647
	estimatedLength += len(s.Seq) * 2 // approximate size for each DNA base
	estimatedLength += len(s.Qual)
	estimatedLength += len(s.Extra)

	return estimatedLength
}
