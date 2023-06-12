package bed

import (
	"io"
)

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

// GetChrom returns the chrom of the bed struct.
func (b Bed) GetChrom() string {
	return b.Chrom
}

// GetChromStart returns the starting coordinates of the bed struct.
func (b Bed) GetChromStart() int {
	return b.ChromStart
}

// GetChromEnd returns the end coordinates of the bed struct.
func (b Bed) GetChromEnd() int {
	return b.ChromEnd
}

// UpdateCoord will return a copy of the bed with the modified coordinates.
func (b Bed) UpdateCoord(c string, start int, end int) interface{} {
	b.Chrom = c
	b.ChromStart = start
	b.ChromEnd = end
	return b
}

// WriteToFileHandle will write a bed struct to the
// io.Writer.
func (b Bed) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, b) //adaptive field writing seems most flexible for the method
}
