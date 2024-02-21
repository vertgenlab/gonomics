package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/giraf"
)

func TestPathToSeq(t *testing.T) {
	genome := &GenomeGraph{
		Nodes: []Node{
			{Seq: dna.StringToBases("ACGT")},
			{Seq: dna.StringToBases("TGCA")},
			{Seq: dna.StringToBases("AAAA")},
		},
	}

	path := giraf.Path{
		Nodes:  []uint32{0, 1, 2},
		TStart: 1,
		TEnd:   7,
	}

	expected := dna.StringToBases("CGTTGCAAAAA")
	result := PathToSeq(path, genome)

	if dna.BasesToString(expected) != dna.BasesToString(result) {
		t.Errorf("PathToSeq returned incorrect result. Expected: %s, Got: %s", dna.BasesToString(expected), dna.BasesToString(result))
	}
}
