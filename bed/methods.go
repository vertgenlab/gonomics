package bed

import (
	"github.com/vertgenlab/gonomics/fileio"
	"io"
)

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (b Bed) GetChrom() string {
	return b.Chrom
}

func (b Bed) GetChromStart() int {
	return int(b.ChromStart)
}

func (b Bed) GetChromEnd() int {
	return int(b.ChromEnd)
}

func (b Bed) UpdateLift(c string, start int, end int) {
	b.Chrom = c
	b.ChromStart = start
	b.ChromEnd = end
}

type BedSlice []*Bed

func (b BedSlice) Len() int { return len(b) }

func (b BedSlice) Swap(i, j int) { b[i], b[j] = b[j], b[i] }

func (b *BedSlice) Push(x interface{}) {
	answer := x.(*Bed)
	*b = append(*b, answer)
}

func (b *BedSlice) Pop() interface{} {
	oldQueue := *b
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*b = oldQueue[:n-1]
	return answer
}

func (b BedSlice) Write(file string) {
	a := make([]Bed, len(b))
	for i := range b {
		a[i] = *b[i]
	}
	Write(file, a)
}

func (b Bed) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, b) //adaptive field writing seems most flexible for the method
}

func (b *Bed) NextRealRecord(file *fileio.EasyReader) bool {
	var done bool
	var next Bed
	for next.Chrom == "" && !done {
		next, done = NextBed(file)
	}
	if done {
		return done
	}
	*b = next
	return done
}

func (b *Bed) Copy() interface{} {
	var answer *Bed = new(Bed)
	*answer = *b
	return answer
}

type ByGenomicCoordinates struct {
	BedSlice
}

func (g ByGenomicCoordinates) Less(i, j int) bool {
	// First sort criteria is chromosome
	if g.BedSlice[i].GetChrom() < g.BedSlice[j].GetChrom() {
		return true
	} else if g.BedSlice[i].GetChrom() == g.BedSlice[j].GetChrom() {
		// If chroms are equal then sort by start position
		if g.BedSlice[i].GetChromStart() < g.BedSlice[j].GetChromStart() {
			return true
		} else if g.BedSlice[i].GetChromStart() == g.BedSlice[j].GetChromStart() {
			// If start positions are equal then the shorter region wins
			if g.BedSlice[i].GetChromEnd() < g.BedSlice[j].GetChromEnd() {
				return true
			}
		}
	}
	return false
}
