package bed

import "github.com/vertgenlab/gonomics/fileio"

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (b *Bed) GetChrom() string {
	return b.Chrom
}

func (b *Bed) GetChromStart() int {
	return int(b.ChromStart)
}

func (b *Bed) GetChromEnd() int {
	return int(b.ChromEnd)
}

type ByGenomicCoordinates []*Bed

func (g ByGenomicCoordinates) Len() int { return len(g) }

func (g ByGenomicCoordinates) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func (g ByGenomicCoordinates) Less(i, j int) bool {
	// First sort criteria is chromosome
	if g[i].GetChrom() < g[j].GetChrom() {
		return true
	} else if g[i].GetChrom() == g[j].GetChrom() {
		// If chroms are equal then sort by start position
		if g[i].GetChromStart() < g[j].GetChromStart() {
			return true
		} else if g[i].GetChromStart() == g[j].GetChromStart() {
			// If start positions are equal then the shorter region wins
			if g[i].GetChromEnd() < g[j].GetChromEnd() {
				return true
			}
		}
	}
	return false
}

func (g *ByGenomicCoordinates) Push(x interface{}) {
	answer := x.(*Bed)
	*g = append(*g, answer)
}

func (g *ByGenomicCoordinates) Pop() interface{} {
	oldQueue := *g
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*g = oldQueue[:n-1]
	return answer
}

func (g ByGenomicCoordinates) Write(file string) {
	Write(file, g, 7)
}

func (b *Bed) WriteToFileHandle(file *fileio.EasyWriter) {
	//TODO: write max fields that are non-nil?
	WriteToFileHandle(file, b, 7)
}

func (b *Bed) NextLine(file *fileio.EasyReader) bool {
	var done bool
	var next *Bed
	next, done = NextBed(file)
	*b = *next
	return done
}