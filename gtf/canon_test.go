package gtf

import (
	"math/rand"
	"testing"
)

func TestCanon(t *testing.T) {
	gtf := Read("testdata/gtfFileTest.gtf")
	for _, i := range gtf {
		SortTranscripts(i)
		if CdsLength(i.Transcripts[0]) != 219 {
			t.Errorf("ERROR: Problem calculating CDS Length")
		}
		if CdnaLength(i.Transcripts[0]) != 1237 {
			t.Errorf("ERROR: Problem calculating cDNA Length")
		}
	}
}

func ShuffleTranscripts(m map[string]*Gene) {
	for _, g := range m {
		rand.Shuffle(len(g.Transcripts), func(i, j int) { g.Transcripts[i], g.Transcripts[j] = g.Transcripts[j], g.Transcripts[i] })
	}
}

func BenchmarkSortTranscripts(b *testing.B) {
	gtf := Read("testdata/long.gtf")
	for i := 0; i < b.N; i++ {
		SortAllTranscripts(gtf)
		ShuffleTranscripts(gtf)
	}
}
