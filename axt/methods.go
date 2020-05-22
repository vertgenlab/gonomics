package axt

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
