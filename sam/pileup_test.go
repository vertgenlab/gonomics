package sam

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"runtime"
	"testing"
	"time"
)

func BenchmarkPileupLinkedList(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		reads, header := GoReadToChan("testdata/peak.bam")
		b.StartTimer()
		GoPileup(reads, header, false, nil, nil)
	}
}

func TestReadBeginsWithInsertion(t *testing.T) {
	reads, header := GoReadToChan("testdata/memtest.bam")
	piles := GoPileup(reads, header, false, nil, nil)
	m := new(runtime.MemStats)
	ticker := time.NewTicker(100 * time.Millisecond)

	var open bool
	for {
		select {
		case <-ticker.C:
			runtime.ReadMemStats(m)
			if m.TotalAlloc > 100000000 { // 100MB
				t.Error("pileup taking too much memory, potential memory leak, see sam/pileup.go func getPile for potential culprit")
				t.FailNow()
			}
		case _, open = <-piles:
			if !open {
				return
			}
		}
	}
}

func TestSyncPileup(t *testing.T) {
	alphaReads, alphaHeader := GoReadToChan("testdata/peak.bam")
	alpha := GoPileup(alphaReads, alphaHeader, false, nil, nil)
	betaReads, betaHeader := GoReadToChan("testdata/peak.bam")
	beta := GoPileup(betaReads, betaHeader, false, nil, nil)
	syncedPilesChan := GoSyncPileups(alpha, beta)
	for syncedPiles := range syncedPilesChan {
		if len(syncedPiles) != 2 {
			t.Error("problem syncing piles")
		}
		if syncedPiles[0].CountF != syncedPiles[1].CountF || syncedPiles[0].CountR != syncedPiles[1].CountR {
			t.Error("problem syncing piles")
		}
	}
}

func TestPeakPileup(t *testing.T) {
	reads, header := GoReadToChan("testdata/peak.bam")
	piles := GoPileup(reads, header, false, nil, nil)
	for pile := range piles {
		switch pile.Pos {
		case 130592024:
			if pile.CountF[dna.A]+pile.CountR[dna.A] != 243 ||
				pile.InsCountF["GAAG"]+pile.InsCountR["GAAG"] != 2 ||
				pile.CountF[dna.Gap]+pile.CountR[dna.Gap] != 4 {
				t.Error("problem with pileup")
			}

		case 130592002:
			if pile.CountF[dna.A]+pile.CountR[dna.A] != 238 {
				t.Error("problem with pileup")
			}

		case 130592001:
			if pile.CountF[dna.G]+pile.CountR[dna.G] != 239 || pile.CountF[dna.C]+pile.CountR[dna.C] != 1 {
				t.Error("problem with pileup")
			}

		case 130592072:
			if pile.CountF[dna.G]+pile.CountR[dna.G] != 237 || pile.CountF[dna.C]+pile.CountR[dna.C] != 1 {
				t.Error("problem with pileup")
			}

		case 130592095:
			if pile.CountF[dna.C]+pile.CountR[dna.C] != 234 {
				t.Error("problem with pileup")
			}
		}
	}
}

func TestRandPileup(t *testing.T) {
	reads, header := GoReadToChan("testdata/rand.bam")
	refmap := chromInfo.SliceToMap(header.Chroms)
	piles := GoPileup(reads, header, false, nil, nil)
	for pile := range piles {
		switch {
		case pile.Pos == 130592072 && pile.RefIdx == refmap["chr9"].Order:
			if pile.CountF[dna.G]+pile.CountR[dna.G] != 2 {
				t.Error("problem with pileup")
			}

		case pile.Pos == 31624960 && pile.RefIdx == refmap["chr18"].Order:
			if pile.CountF[dna.G]+pile.CountR[dna.G] != 2 {
				t.Error("problem with pileup")
			}

		case pile.Pos == 24954689 && pile.RefIdx == refmap["chrX"].Order:
			if pile.CountF[dna.C]+pile.CountR[dna.C] != 2 {
				t.Error("problem with pileup")
			}

		case pile.Pos == 45795462 && pile.RefIdx == refmap["chr12"].Order:
			if pile.CountF[dna.T]+pile.CountR[dna.T] != 1 {
				t.Error("problem with pileup")
			}

		case pile.Pos == 91864875 && pile.RefIdx == refmap["chr7"].Order:
			if pile.CountF[dna.T]+pile.CountR[dna.T] != 2 {
				t.Error("problem with pileup")
			}
		}
	}
}

var p1 = Sam{
	RName: "ref",
	Pos:   1,
	Cigar: cigar.FromString("2M1I1M1S"),
	Seq:   dna.StringToBases("TTAGA"),
}

var p2 = Sam{
	RName: "ref",
	Pos:   2,
	Cigar: cigar.FromString("5H2M5D3M"), // 5bp deletion at position 3
	Seq:   dna.StringToBases("AAAAG"),
}

var p3 = Sam{
	RName: "ref",
	Pos:   5,
	Cigar: cigar.FromString("10M"),
	Seq:   dna.StringToBases("CCCCCCCCCC"),
}

var p4 = Sam{
	RName: "ref",
	Pos:   14,
	Cigar: cigar.FromString("5M"),
	Seq:   dna.StringToBases("GGGGG"),
}

var p5 = Sam{
	RName: "ref",
	Pos:   20,
	Flag:  81, // paired + reverse strand + first in pair
	Cigar: cigar.FromString("5M"),
	Seq:   dna.StringToBases("GCGCG"),
}

var p6 = Sam{
	RName: "ref",
	Pos:   20,
	Flag:  161, // paired + mate reverse strand + second in pair
	Cigar: cigar.FromString("5I5M5I"),
	Seq:   dna.StringToBases("AAAAAGCCCGAAAAA"),
}

var pHeader = Header{
	Chroms: []chromInfo.ChromInfo{
		{Name: "ref", Order: 0, Size: 30},
	},
	Metadata: Metadata{
		SortOrder: []SortOrder{Coordinate},
	},
}

func TestMultiBaseDeletion(t *testing.T) {
	reads := make(chan Sam, 6)
	reads <- p1
	reads <- p2
	reads <- p3
	reads <- p4
	reads <- p5
	reads <- p6
	close(reads)

	piles := GoPileup(reads, pHeader, true, nil, nil)

	actual := make(map[uint32]Pile) // keyed by position
	for p := range piles {
		actual[p.Pos] = copyPile(p)
	}

	if len(actual) != 30 { // len of chrom, includeNoData is true
		t.Error("problem with includeNoData")
	}

	if len(actual[4].DelCountF) != 1 || len(actual[4].DelCountR) != 0 {
		t.Error("problem with deletions")
	}

	if actual[4].DelCountF[5] != 1 {
		t.Error("problem with deletions")
	}

	var i uint32
	for i = 0; i < 20; i++ {
		if actual[i].CountR != [13]int{} { // should be empty until 20
			t.Error("problem with read strand")
		}
	}

	if actual[20].CountF != actual[20].CountR {
		t.Error("problem with read strand")
	}
}

func copyPile(p Pile) Pile {
	var n Pile
	n = p

	if p.InsCountF != nil {
		n.InsCountF = make(map[string]int)
		for key, val := range p.InsCountF {
			n.InsCountF[key] = val
		}
	}

	if p.InsCountR != nil {
		n.InsCountR = make(map[string]int)
		for key, val := range p.InsCountR {
			n.InsCountR[key] = val
		}
	}

	if p.DelCountF != nil {
		n.DelCountF = make(map[int]int)
		for key, val := range p.DelCountF {
			n.DelCountF[key] = val
		}
	}

	if p.DelCountR != nil {
		n.DelCountR = make(map[int]int)
		for key, val := range p.DelCountR {
			n.DelCountR[key] = val
		}
	}
	return n
}

var PileConsensusTests = []struct {
	p                  Pile
	c                  Consensus
	substitutionsOnly  bool
	insertionThreshold float64
}{
	{p: Pile{
		CountF:    [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{"AAT": 2, "AC": 1},
		InsCountR: map[string]int{"AAT": 1, "AC": 3},
		DelCountF: map[int]int{6: 16, 5: 2},
		DelCountR: map[int]int{6: 19, 5: 1},
	},
		c: Consensus{
			Base:      0,
			Insertion: []dna.Base{},
			Deletion:  6,
			Type:      Deletion,
		},
		substitutionsOnly:  false,
		insertionThreshold: 0.1,
	},
	{p: Pile{
		CountF:    [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{},
		InsCountR: map[string]int{},
		DelCountF: map[int]int{},
		DelCountR: map[int]int{},
	},
		c: Consensus{
			Base:      0,
			Insertion: []dna.Base{},
			Deletion:  0,
			Type:      Undefined,
		},
		substitutionsOnly:  false,
		insertionThreshold: 0.1,
	},
	{p: Pile{
		CountF:    [13]int{0, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 17, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{"AAT": 2, "AC": 1},
		InsCountR: map[string]int{"AAT": 1, "AC": 3},
		DelCountF: map[int]int{6: 2, 5: 2},
		DelCountR: map[int]int{6: 4, 5: 1},
	},
		c: Consensus{
			Base:      1,
			Insertion: []dna.Base{},
			Deletion:  0,
			Type:      Base,
		},
		substitutionsOnly:  false,
		insertionThreshold: 0.1,
	},
	{p: Pile{
		CountF:    [13]int{0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{"AACTT": 42, "AC": 1},
		InsCountR: map[string]int{"AACTT": 49, "AC": 3},
		DelCountF: map[int]int{6: 2, 5: 2},
		DelCountR: map[int]int{6: 4, 5: 1},
	},
		c: Consensus{
			Base:      1,
			Insertion: dna.StringToBases("AACTT"),
			Deletion:  0,
			Type:      Insertion,
		},
		substitutionsOnly:  false,
		insertionThreshold: 0.1,
	},
	{p: Pile{
		CountF:    [13]int{0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{"AACTT": 42, "AC": 1},
		InsCountR: map[string]int{"AACTT": 49, "AC": 3},
		DelCountF: map[int]int{6: 2, 5: 2},
		DelCountR: map[int]int{6: 4, 5: 1},
	},
		c: Consensus{
			Base:      1,
			Insertion: []dna.Base{},
			Deletion:  0,
			Type:      Base,
		},
		substitutionsOnly:  true,
		insertionThreshold: 0.1,
	},
}

func TestPileConsensus(t *testing.T) {
	var answer Consensus
	for _, v := range PileConsensusTests {
		answer = PileConsensus(v.p, v.substitutionsOnly, v.insertionThreshold)
		if answer.Type != v.c.Type {
			t.Errorf("Error in PileConsensus. Observed consensus type was not as expected.")
		}
		switch v.c.Type {
		case Base:
			if v.c.Base != answer.Base {
				t.Errorf("Error in PileConsensus. Observed consensus base did not match expected.")
			}
		case Insertion:
			if dna.BasesToString(v.c.Insertion) != dna.BasesToString(answer.Insertion) {
				t.Errorf("Error in PileConsensus. Observed consensus insertion did not match expected.")
			}
		case Deletion:
			if v.c.Deletion != answer.Deletion {
				t.Errorf("Error in PileConsensus. Observed consensus deletion did not match expected.")
			}
			// we already checked type so no need to check anything else for undefined consensus.
		}
	}
}
