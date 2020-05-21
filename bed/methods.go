package bed

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

type byGenomicCoordinates []*Bed

func (g byGenomicCoordinates) Len() int { return len(g) }

func (g byGenomicCoordinates) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func (g byGenomicCoordinates) Less(i, j int) bool {
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

func (g *byGenomicCoordinates) Push(x interface{}) {
	answer := x.(*Bed)
	*g = append(*g, answer)
}

func (g *byGenomicCoordinates) Pop() interface{} {
	oldQueue := *g
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*g = oldQueue[:n-1]
	return answer
}