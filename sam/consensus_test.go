package sam

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

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
		InsCountR: map[string]int{"AAT": 1, "AC": 1},
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
		CountF:    [13]int{0, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		CountR:    [13]int{0, 17, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		InsCountF: map[string]int{"AAT": 2, "AC": 1},
		InsCountR: map[string]int{"AAT": 1, "AC": 4},
		DelCountF: map[int]int{6: 2, 5: 2},
		DelCountR: map[int]int{6: 4, 5: 1},
	},
		c: Consensus{
			Base:      1,
			Insertion: dna.StringToBases("AC"),
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
			t.Errorf("Error in PileConsensus. Observed consensus type: %v. was not as expected: %v.", answer.Type, v.c.Type)
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
