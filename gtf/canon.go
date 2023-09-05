package gtf

import "sort"

// CdnaLength returns the length of the cDNA in nucleotides.
func CdnaLength(t *Transcript) int {
	var answer int = 0
	for i := 0; i < len(t.Exons); i++ {
		answer += t.Exons[i].End - t.Exons[i].Start + 1
	}
	return answer
}

// CdsLength returns the length of the Cds in nucleotides (after splicing).
func CdsLength(t *Transcript) int {
	var answer int = 0
	for i := 0; i < len(t.Exons); i++ {
		if t.Exons[i].Cds != nil {
			answer += t.Exons[i].Cds.End - t.Exons[i].Cds.Start + 1
		}
	}
	return answer
}

func isLonger(i, j *Transcript) bool {
	iLen := CdsLength(i)
	jLen := CdsLength(j)
	return iLen > jLen || (iLen == jLen && CdnaLength(i) > CdnaLength(j))
}

// SortTranscripts sorts the longest transcript to the front so that the canonical/longest transcript is always g.Transcripts[0].
func SortTranscripts(g *Gene) {
	sort.Slice(g.Transcripts, func(i, j int) bool { return isLonger(g.Transcripts[i], g.Transcripts[j]) })
}

// SortAllTranscripts applies SortTranscripts to every value in the map
func SortAllTranscripts(m map[string]*Gene) {
	for _, g := range m {
		SortTranscripts(g)
	}
}

// MoveCanonicalToZero does a single iteration of bubble sort to move the longest/canonical transcript to the first position in the slice.
// This is faster than SortTranscripts.
func MoveCanonicalToZero(g *Gene) {
	for i := 1; i < len(g.Transcripts); i++ {
		if isLonger(g.Transcripts[i], g.Transcripts[0]) {
			g.Transcripts[0], g.Transcripts[i] = g.Transcripts[i], g.Transcripts[0]
		}
	}
}

// MoveAllCanonicalToZero applies MoveCanonicalToZero to every value in the map
func MoveAllCanonicalToZero(m map[string]*Gene) {
	for _, g := range m {
		MoveCanonicalToZero(g)
	}
}
