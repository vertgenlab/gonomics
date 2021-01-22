package dna

import "testing"

var alpha []Base = StringToBases("ACTGacgtAAACC--ACacgnnnNNNactg")
var beta []Base = StringToBases("ACGGacgtAATCC--ACacgnnnNCNaCtg")

func TestCount(t *testing.T) {
	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount := Count(alpha)
	if ACount != 5 || CCount != 4 || GCount != 1 || TCount != 1 || NCount != 3 || aCount != 3 || cCount != 3 ||
		gCount != 3 || tCount != 2 || nCount != 3 || gapCount != 2 {
		t.Errorf("problem with Count")
	}
}

func TestCountMask(t *testing.T) {
	umC, mC, gC := CountMask(alpha)
	if umC != 14 || mC != 14 || gC != 2 {
		t.Errorf("problem with CountMask")
	}
}

func TestCountGaps(t *testing.T) {
	if CountGaps(alpha) != 2 {
		t.Errorf("problem with CountGaps")
	}
}

func TestDist(t *testing.T) {
	actual, err := Dist(alpha, beta)
	if actual != 4 || err != nil {
		t.Errorf("problem with Dist")
	}
}

func TestCountBase(t *testing.T) {
	actual := CountBase(alpha, N)
	if actual != 3 {
		t.Errorf("problem with CountBase")
	}
}

func TestDefineBase(t *testing.T) {
	if !DefineBase(A) || !DefineBase(C) || !DefineBase(G) || !DefineBase(T) ||
		!DefineBase(LowerA) || !DefineBase(LowerC) || !DefineBase(LowerG) ||
		!DefineBase(LowerT) || DefineBase(N) || DefineBase(LowerN) ||
		DefineBase(Gap) || DefineBase(Dot) || DefineBase(Nil) {
		t.Errorf("problem with DefineBase")
	}
}