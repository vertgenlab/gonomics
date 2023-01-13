package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

type seedTest struct {
	SeqA             []dna.Base
	SmallA           *ThreeBit
	SeqB             []dna.Base
	SmallB           *ThreeBit
	StartA           int
	StartB           int
	TrueMatchesLeft  int
	TrueMatchesRight int
}

var longDnaOne []dna.Base = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGG")
var longDnaTwo []dna.Base = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var longDnaThree []dna.Base = dna.StringToBases("CTGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var longDnaFour []dna.Base = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")

var shortDnaOne []dna.Base = dna.StringToBases("CCCCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGG")
var shortDnaTwo []dna.Base = dna.StringToBases("ACCTACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var shortDnaThree []dna.Base = dna.StringToBases("T")
var shortDnaFour []dna.Base = dna.StringToBases("T")

var seedTestsLong []seedTest = []seedTest{{SeqA: longDnaOne,
	SmallA:           NewThreeBit(longDnaOne, PaddingOne),
	SeqB:             longDnaTwo,
	SmallB:           NewThreeBit(longDnaTwo, PaddingTwo),
	StartA:           100,
	StartB:           100,
	TrueMatchesLeft:  101,
	TrueMatchesRight: 43,
},
	{SeqA: longDnaThree,
		SmallA:           NewThreeBit(longDnaThree, PaddingOne),
		SeqB:             longDnaFour,
		SmallB:           NewThreeBit(longDnaFour, PaddingTwo),
		StartA:           90,
		StartB:           90,
		TrueMatchesLeft:  90,
		TrueMatchesRight: 54,
	},
}

var seedTestsShort []seedTest = []seedTest{{SeqA: shortDnaOne,
	SmallA:           NewThreeBit(shortDnaOne, PaddingOne),
	SeqB:             shortDnaTwo,
	SmallB:           NewThreeBit(shortDnaTwo, PaddingTwo),
	StartA:           1,
	StartB:           1,
	TrueMatchesLeft:  1,
	TrueMatchesRight: 2,
},
	{SeqA: shortDnaThree,
		SmallA:           NewThreeBit(shortDnaThree, PaddingOne),
		SeqB:             shortDnaFour,
		SmallB:           NewThreeBit(shortDnaFour, PaddingTwo),
		StartA:           0,
		StartB:           0,
		TrueMatchesLeft:  1,
		TrueMatchesRight: 1,
	},
}

func currentMethodRight(seqA []dna.Base, startA int, seqB []dna.Base, startB int) int {
	var matches int = 0
	for i, j := startA, startB; i < len(seqA) && j < len(seqB) && seqA[i] == seqB[j]; i, j = i+1, j+1 {
		matches++
	}
	return matches
}

func currentMethodLeft(seqA []dna.Base, startA int, seqB []dna.Base, startB int) int {
	var matches int = 0
	for i, j := startA, startB; i >= 0 && j >= 0 && seqA[i] == seqB[j]; i, j = i-1, j-1 {
		matches++
	}
	return matches
}

func TestCounting(t *testing.T) {
	var matches int
	seedTests := append(seedTestsShort, seedTestsLong...)
	for _, data := range seedTests {
		matches = countLeftMatches(data.SmallA.Seq, data.StartA, data.SmallB.Seq, data.StartB)
		if matches != data.TrueMatchesLeft {
			t.Errorf("Error: found %d left matches, but was expecting %d\n", matches, data.TrueMatchesLeft)
		}
		matches = countRightMatches(data.SmallA.Seq, data.StartA, data.SmallA.Len, data.SmallB.Seq, data.StartB, data.SmallB.Len)
		if matches != data.TrueMatchesRight {
			t.Errorf("Error: found %d right matches, but was expecting %d\n", matches, data.TrueMatchesRight)
		}
	}
}

func BenchmarkCurrentShortLeft(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsShort {
			matches = currentMethodLeft(data.SeqA, data.StartA, data.SeqB, data.StartB)
			if matches != data.TrueMatchesLeft {
				b.Errorf("Error: wrong number of matches.\n")
			}
		}
	}
}

func BenchmarkCurrentLongLeft(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsLong {
			matches = currentMethodLeft(data.SeqA, data.StartA, data.SeqB, data.StartB)
			if matches != data.TrueMatchesLeft {
				b.Errorf("Error: wrong number of left matches.  Expected=%d Calculated=%d\n", data.TrueMatchesLeft, matches)
			}
		}
	}
}

func BenchmarkCurrentShortRight(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsShort {
			matches = currentMethodRight(data.SeqA, data.StartA, data.SeqB, data.StartB)
			if matches != data.TrueMatchesRight {
				b.Errorf("Error: wrong number of right matches.  Expected=%d Calculated=%d\n", data.TrueMatchesRight, matches)
			}
		}
	}
}

func BenchmarkCurrentLongRight(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsLong {
			matches = currentMethodRight(data.SeqA, data.StartA, data.SeqB, data.StartB)
			if matches != data.TrueMatchesRight {
				b.Errorf("Error: wrong number of matches.\n")
			}
		}
	}
}

func BenchmarkNewShortLeft(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsShort {
			matches = countLeftMatches(data.SmallA.Seq, data.StartA, data.SmallB.Seq, data.StartB)
			if matches != data.TrueMatchesLeft {
				b.Errorf("Error: wrong number of matches.\n")
			}
		}
	}
}

func BenchmarkNewLongLeft(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsLong {
			matches = countLeftMatches(data.SmallA.Seq, data.StartA, data.SmallB.Seq, data.StartB)
			if matches != data.TrueMatchesLeft {
				b.Errorf("Error: wrong number of matches.\n")
			}
		}
	}
}

func BenchmarkNewShortRight(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsShort {
			matches = countRightMatches(data.SmallA.Seq, data.StartA, data.SmallA.Len, data.SmallB.Seq, data.StartB, data.SmallB.Len)
			if matches != data.TrueMatchesRight {
				b.Errorf("Error: wrong number of matches\n")
			}
		}
	}
}

func BenchmarkNewLongRight(b *testing.B) {
	var matches int

	for n := 0; n < b.N; n++ {
		for _, data := range seedTestsLong {
			matches = countRightMatches(data.SmallA.Seq, data.StartA, data.SmallA.Len, data.SmallB.Seq, data.StartB, data.SmallB.Len)
			if matches != data.TrueMatchesRight {
				b.Errorf("Error: wrong number of matches\n")
			}
		}
	}
}
