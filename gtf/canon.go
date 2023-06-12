package gtf

import "sort"

// Returns length of cDNA in nucleotides.
func CdnaLength(t *Transcript) int {
	var answer int = 0
	for i := 0; i < len(t.Exons); i++ {
		answer += t.Exons[i].End - t.Exons[i].Start + 1
	}
	return answer
}

// Returns length of Cds in nucleotides.
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

// Sorts so that the canonical transcript is always g.Transcripts[0].
func SortTranscripts(g *Gene) {
	sort.Slice(g.Transcripts, func(i, j int) bool { return isLonger(g.Transcripts[i], g.Transcripts[j]) })
}

func SortAllTranscripts(m map[string]*Gene) {
	for _, g := range m {
		SortTranscripts(g)
	}
}

// Moves canonical to zero without sorting all transcripts. Faster than SortTranscripts.
func MoveCanonicalToZero(g *Gene) {
	for i := 1; i < len(g.Transcripts); i++ {
		if isLonger(g.Transcripts[i], g.Transcripts[0]) {
			g.Transcripts[0], g.Transcripts[i] = g.Transcripts[i], g.Transcripts[0]
		}
	}
}

func MoveAllCanonicalToZero(m map[string]*Gene) {
	for _, g := range m {
		MoveCanonicalToZero(g)
	}
}
