package bed

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
)

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

// GetChrom returns the chrom of the bed struct
func (b Bed) GetChrom() string {
	return b.Chrom
}

// GetChromStart returns the starting coordinates of the bed struct
func (b Bed) GetChromStart() int {
	return b.ChromStart
}

// GetChromEnd returns the end coordinates of the bed struct
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

type BedSlice []*Bed

// Len returns total bed entries within the input bed slice
func (b BedSlice) Len() int { return len(b) }

// Swap will switch the values of two bed structs inside
// the slice of beds.
func (b BedSlice) Swap(i, j int) { b[i], b[j] = b[j], b[i] }

// Push and pop() satisfy the interface for the heap method.
// Push pushes to the heap, adding the element to the top of the heap.
func (b *BedSlice) Push(x interface{}) {
	answer := x.(*Bed)
	*b = append(*b, answer)
}

// Pop and push() satisfy the interface for the heap method.
// Pop returns the value at the top to the heap, in the process removing the
// element from the heap.
func (b *BedSlice) Pop() interface{} {
	oldQueue := *b
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*b = oldQueue[:n-1]
	return answer
}

// Write will take a bed slice of bed structs and
// write it to an input file name type string
func (b BedSlice) Write(file string) {
	var err error
	f := fileio.EasyCreate(file)
	for i := range b {
		WriteBed(f, *b[i])
	}
	err = f.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle will write a bed struct to the
// io.Writer
func (b Bed) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, b) //adaptive field writing seems most flexible for the method
}

// NextRealRecord keeps track of if we are done reading through
// the bed by returning a bool and sets b equal to the next bed entry in the file.
func (b *Bed) NextRealRecord(file *fileio.EasyReader) bool {
	var done bool
	var next Bed
	for next.FieldsInitialized == 0 && !done {
		next, done = NextBed(file)
	}
	if done {
		return done
	}
	*b = next
	return done
}

// Copy will return a copy of the bed
func (b *Bed) Copy() interface{} {
	var answer *Bed = new(Bed)
	*answer = *b
	return answer
}

type ByGenomicCoordinates struct {
	BedSlice
}

// Less sorts first by chromosome, then by genomic coordinates starting with the start coordinate
// and breaking a tie with the end coordinates. Less satisfies the sort.Slice
// interface.
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
