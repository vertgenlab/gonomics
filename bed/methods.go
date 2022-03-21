package bed

import (
	"io"
)

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (b Bed) GetChrom() string {
	return b.Chrom
}

func (b Bed) GetChromStart() int {
	return b.ChromStart
}

func (b Bed) GetChromEnd() int {
	return b.ChromEnd
}

func (b Bed) UpdateCoord(c string, start int, end int) interface{} {
	b.Chrom = c
	b.ChromStart = start
	b.ChromEnd = end
	return b
}

func (b Bed) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, b) //adaptive field writing seems most flexible for the method
}
