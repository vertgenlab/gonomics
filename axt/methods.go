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

type AxtSlice []*Axt

func (a AxtSlice) Len() int { return len(a) }

func (a AxtSlice) Swap(i, j int) { a[i], a[j] = a[j], a[i] }

func (a *AxtSlice) Push(x interface{}) {
	answer := x.(*Axt)
	*a = append(*a, answer)
}

func (a *AxtSlice) Pop() interface{} {
	oldQueue := *a
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*a = oldQueue[:n-1]
	return answer
}

func (a AxtSlice) Write(file string) {
	Write(file, a)
}

func (a *Axt) WriteToFileHandle(file *fileio.EasyWriter) {
	//TODO: what to do with alnNumber???
	WriteToFileHandle(file, a, 0)
}

func (a *Axt) NextRealRecord(file *fileio.EasyReader) bool {
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

type ByGenomicCoordinates struct {
	AxtSlice
}

func (g ByGenomicCoordinates) Less(i, j int) bool {
	// First sort criteria is chromosome
	if g.AxtSlice[i].GetChrom() < g.AxtSlice[j].GetChrom() {
		return true
	} else if g.AxtSlice[i].GetChrom() == g.AxtSlice[j].GetChrom() {
		// If chroms are equal then sort by start position
		if g.AxtSlice[i].GetChromStart() < g.AxtSlice[j].GetChromStart() {
			return true
		} else if g.AxtSlice[i].GetChromStart() == g.AxtSlice[j].GetChromStart() {
			// If start positions are equal then the shorter region wins
			if g.AxtSlice[i].GetChromEnd() < g.AxtSlice[j].GetChromEnd() {
				return true
			}
		}
	}
	return false
}