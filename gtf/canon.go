package gtf

import "sort"

// Returns length of cDNA in nucleotides
func CdnaLength(t *Transcript) int {
	var answer int = 0
	for i := 0; i < len(t.Exons); i++ {
		answer += t.Exons[i].End - t.Exons[i].Start
	}
	return answer
}

// Returns length of CDS in nucleotides
func CdsLength(t *Transcript) int {
	var answer int = 0
	for i := 0; i < len(t.Exons); i++ {
		if t.Exons[i].Cds != nil {
			answer += t.Exons[i].Cds.End - t.Exons[i].Cds.Start
		}
	}
	return answer
}

// Sorts so that the cannonical transcript is always g.Transcripts[0]
func SortTranscripts(g *Gene) {
	var dif int

	less := func(i, j int) bool {
		dif = CdsLength(g.Transcripts[i]) - CdsLength(g.Transcripts[j])
		switch {
		case dif < 0:
			return false
		case dif > 0:
			return true
		default:
			return CdnaLength(g.Transcripts[i]) > CdnaLength(g.Transcripts[j])
		}
	}

	sort.Slice(g.Transcripts, less)
}

func SortAllTranscripts(m map[string]*Gene) {
	for _, g := range m {
		SortTranscripts(g)
	}
}
