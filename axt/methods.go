package axt

import "github.com/vertgenlab/gonomics/fileio"

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (a *Axt) GetChrom() string {
	return a.RName
}

// to conform with bed standards the startpos will be zero base and the endpos will be 1 base

func (a *Axt) GetChromStart() int {
	return int(a.RStart - 1)
}

func (a *Axt) GetChromEnd() int {
	return int(a.REnd)
}

type ByGenomicCoordinates []*Axt

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
	answer := x.(*Axt)
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
	Write(file, g)
}

func (a *Axt) WriteToFileHandle(file *fileio.EasyWriter) {
	//TODO: what to do with alnNumber???
	WriteToFileHandle(file, a, 0)
}

func (a *Axt) NextRealLine(file *fileio.EasyReader) bool {
	var done bool
	var next *Axt
	for next == nil && !done {
		next, done = NextAxt(file)
	}
	if done {
		return done
	}
	*a = *next
	return done
}

func (a *Axt) Copy(to *interface{}) {
	var answer *Axt = new(Axt)
	*answer = *a
	*to = answer
}
