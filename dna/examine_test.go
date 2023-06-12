package dna

import "testing"

var alpha = StringToBases("ACTGacgtAAACC--ACacgnnnNNNactg")
var beta = StringToBases("ACGGacgtAATCC--ACacgnnnNCNaCtg")

func TestCount(t *testing.T) {
	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount := Count(alpha)
	if ACount != 5 || CCount != 4 || GCount != 1 || TCount != 1 || NCount != 3 || aCount != 3 || cCount != 3 ||
		gCount != 3 || tCount != 2 || nCount != 3 || gapCount != 2 {
		t.Errorf("problem with Count")
	}

	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount = Count(beta)
	if ACount != 4 || CCount != 6 || GCount != 2 || TCount != 1 || NCount != 2 ||
		aCount != 3 || cCount != 2 || gCount != 3 || tCount != 2 || nCount != 3 ||
		gapCount != 2 {
		t.Errorf("problem with Count")
	}
}

func TestCountMask(t *testing.T) {
	umC, mC, gC := CountMask(alpha)
	if umC != 14 || mC != 14 || gC != 2 {
		t.Error("problem with CountMask(alpha)")
	}

	umC, mC, gC = CountMask(beta)
	if umC != 15 || mC != 13 || gC != 2 {
		t.Error("problem with CountMask(beta)")
	}
}

func TestCountGaps(t *testing.T) {
	if CountGaps(alpha) != 2 {
		t.Errorf("problem with CountGaps")
	}
}

func TestDist(t *testing.T) {
	actual := Dist(alpha, beta)
	if actual != 4 {
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
