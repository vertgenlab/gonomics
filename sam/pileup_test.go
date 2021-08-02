package sam

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var p1 = Sam{
	RName: "ref",
	Pos:   1,
	Cigar: cigar.FromString("2M1I1M1S"),
	Seq:   dna.StringToBases("TTAGA"),
}

var p2 = Sam{
	RName: "ref",
	Pos:   2,
	Cigar: cigar.FromString("5H2M1D3M"),
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
	Cigar: cigar.FromString("5M"),
	Seq:   dna.StringToBases("GGGGG"),
}

var p6 = Sam{
	RName: "ref",
	Pos:   23,
	Cigar: cigar.FromString("5M"),
	Seq:   dna.StringToBases("GGGGG"),
}

var pRef2 = Sam{
	RName: "ref2",
	Pos:   3,
	Cigar: cigar.FromString("3M"),
	Seq:   dna.StringToBases("GGG"),
}

var pRefMap = map[string]chromInfo.ChromInfo {
	"ref": {Order: 0, Size: 30},
	"ref2": {Order: 1, Size: 10},
}

func TestPileBuffer(t *testing.T) {
	testChan := make(chan Pile, 100)
	pb := newPileBuffer(10)
	updatePile(pb, p1, pRefMap) // intialize new buffer
	updatePile(pb, p2, pRefMap)
	if len(pb.s) != 10 {
		t.Error("problem with pile buffer")
	}
	updatePile(pb, p3, pRefMap) // should expand cap to 20 to make space
	if len(pb.s) != 20 {
		t.Error("problem with pile buffer")
	}
	updatePile(pb, p4, pRefMap) // should fill existing space without cap change
	if len(pb.s) != 20 {
		t.Error("problem with pile buffer")
	}

	pb.sendPassed(p3, false, pRefMap, testChan) // should send pos 1:5 and set new pb.idx
	if len(testChan) != 4 || pb.s[pb.idx].Pos != 5 {
		t.Error("problem with pile buffer sendPassed")
	}

	updatePile(pb, p5, pRefMap) // should wrap to start of buffer
	if pb.s[0].Pos != 21 {
		t.Error("problem with pile buffer wrapping")
	}

	if pb.getIdx(0, 5).Pos != 5 {
		t.Error("did not get right pos")
	}
	if pb.getIdx(0, 19).Pos != 19 {
		t.Error("did not get right pos")
	}
	if pb.getIdx(0, 23).Pos != 23 {
		t.Error("did not get right pos")
	}

	updatePile(pb, p6, pRefMap) // should expand buffer and reorder starting at 0
	if pb.s[0].Pos != 5 || len(pb.s) != 40 {
		t.Error("problem with pile buffer expansion")
	}
}

func TestPileBufferSending(t *testing.T) {
	testChan := make(chan Pile, 100)
	pb := newPileBuffer(10)
	updatePile(pb, p4, pRefMap)
	updatePile(pb, p5, pRefMap) // leave untouched base at pos 19

	pb.sendPassed(pRef2, false, pRefMap, testChan) // should send 10 piles
	if len(testChan) != 10 {
		t.Error("problem with pile buffer sendPassed")
	}

	// include no data tests
	testChan = make(chan Pile, 100)
	pb = newPileBuffer(10)
	updatePile(pb, p4, pRefMap)
	updatePile(pb, p5, pRefMap) // leave untouched base at pos 19

	pb.sendPassed(pRef2, true, pRefMap, testChan) // should send 10 piles

	close(testChan)
	var foundMissing bool
	for i := range testChan {
		if i.Pos == 19 {
			foundMissing = true
		}
	}
	if !foundMissing {
		t.Error("did not find piles for no data")
	}
}