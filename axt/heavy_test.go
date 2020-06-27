package axt

import (
	"testing"
)

func TestAxtHeavy(t *testing.T) {
	// Making sure methods for AxtHeavy overrule methods for Axt
	a := &Axt{
		RName: "R",
		RStart: 0,
		REnd: 10,
		QName: "Q",
		QStart: 100,
		QEnd: 110,
		QStrandPos: false,
	}

	b := AxtHeavy {
		Axt: a,
		Chrom: "Heavy",
		ChromSize: 200,
		ChromStart: 50,
		ChromEnd: 60,
	}

	if b.GetChrom() == a.GetChrom() ||
		b.GetChromStart() == a.GetChromStart() ||
		b.GetChromEnd() == a.GetChromEnd() {
		t.Errorf("ERROR: Methods for Axt overrule methods for AxtHeavy")
	}
}
