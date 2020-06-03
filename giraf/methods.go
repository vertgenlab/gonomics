package giraf

import (
	"github.com/vertgenlab/gonomics/fileio"
)

func (g *Giraf) GetChrom() string {
	return string(g.Path.Nodes[0])
}

func (g *Giraf) GetChromStart() int {
	return g.Path.TStart
}

func (g *Giraf) GetChromEnd() int {
	// BEWARE this returns the end position on the final node
	// There may be a gotcha if you assume GetChromStart() and End()
	// refer to positions on the same chromosome
	// Not sure how best to handle this???
	return g.Path.TEnd
}

type GirafSlice []*Giraf

func (g GirafSlice) Len() int { return len(g) }

func (g GirafSlice) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func (g *GirafSlice) Push(x interface{}) {
	answer := x.(*Giraf)
	*g = append(*g, answer)
}

func (g *GirafSlice) Pop() interface{} {
	oldQueue := *g
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*g = oldQueue[:n-1]
	return answer
}

func (g GirafSlice) Write(filename string) {
	Write(filename, g)
}

func (g *Giraf) WriteToFileHandle(file *fileio.EasyWriter) {
	//TODO: change WriteGirafToFileHandle to operate on fileio.EasyReader instead of os.File
	WriteGirafToFileHandle(file.File, g)
}

func (g *Giraf) NextRealRecord(file *fileio.EasyReader) bool {
	var done bool
	var next *Giraf
	for next == nil && !done {
		next, done = NextGiraf(file)
	}
	if done {
		return true
	}
	*g = *next
	return done
}

func (g *Giraf) Copy() interface{} {
	var answer *Giraf = new(Giraf)
	*answer = *g
	return answer
}

type ByTopologicalNodeOrder struct {
	GirafSlice
}

func (g ByTopologicalNodeOrder) Less(i, j int) bool {
	// First sort criteria is node
	if g.GirafSlice[i].GetChrom() < g.GirafSlice[j].GetChrom() {
		return true
	} else if g.GirafSlice[i].GetChrom() == g.GirafSlice[j].GetChrom() {
		// If start nodes are equal then sort by start position
		if g.GirafSlice[i].GetChromStart() < g.GirafSlice[j].GetChromStart() {
			return true
		} else if g.GirafSlice[i].GetChromStart() == g.GirafSlice[j].GetChromStart() {
			// If start positions are equal then loop through nodes and see if one has priority
			for k := 0; k < len(g.GirafSlice[i].Path.Nodes); k++ {
				if g.GirafSlice[i].Path.Nodes[k] < g.GirafSlice[j].Path.Nodes[k] {
					return true
				}
			}
			// If all nodes match, sort based on longest path
			if len(g.GirafSlice[i].Path.Nodes) < len(g.GirafSlice[j].Path.Nodes) {
				return true
			} else if len(g.GirafSlice[i].Path.Nodes) == len(g.GirafSlice[j].Path.Nodes) {
				// If nodes are equal length, then sort based on the ending position
				if g.GirafSlice[i].GetChromEnd() < g.GirafSlice[j].GetChromEnd() {
					return true
				}
			}
		}
	}
	return false
}
