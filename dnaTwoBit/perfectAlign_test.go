package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

type seedTest struct {
	SeqA             []dna.Base
	SmallA           *TwoBit
	SeqB             []dna.Base
	SmallB           *TwoBit
	StartA           int
	StartB           int
	TrueMatchesLeft  int
	TrueMatchesRight int
}

var longDnaOne, _ = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGG")
var longDnaTwo, _ = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var longDnaThree, _ = dna.StringToBases("CTGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var longDnaFour, _ = dna.StringToBases("ATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")

var shortDnaOne, _ = dna.StringToBases("CCCCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGG")
var shortDnaTwo, _ = dna.StringToBases("ACCTACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGTATGCACTAGTCATACAGT")
var shortDnaThree, _ = dna.StringToBases("T")
var shortDnaFour, _ = dna.StringToBases("T")

var seedTestsLong []seedTest = []seedTest{seedTest{SeqA: longDnaOne,
	SmallA:           NewTwoBit(longDnaOne),
	SeqB:             longDnaTwo,
	SmallB:           NewTwoBit(longDnaTwo),
	StartA:           100,
	StartB:           100,
	TrueMatchesLeft:  101,
	TrueMatchesRight: 43,
},
	seedTest{SeqA: longDnaThree,
		SmallA:           NewTwoBit(longDnaThree),
		SeqB:             longDnaFour,
		SmallB:           NewTwoBit(longDnaFour),
		StartA:           90,
		StartB:           90,
		TrueMatchesLeft:  90,
		TrueMatchesRight: 54,
	},
}

var seedTestsShort []seedTest = []seedTest{seedTest{SeqA: shortDnaOne,
	SmallA:           NewTwoBit(shortDnaOne),
	SeqB:             shortDnaTwo,
	SmallB:           NewTwoBit(shortDnaTwo),
	StartA:           1,
	StartB:           1,
	TrueMatchesLeft:  1,
	TrueMatchesRight: 2,
},
	seedTest{SeqA: shortDnaThree,
		SmallA:           NewTwoBit(shortDnaThree),
		SeqB:             shortDnaFour,
		SmallB:           NewTwoBit(shortDnaFour),
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
		matches = CountLeftMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
		if matches != data.TrueMatchesLeft {
			t.Errorf("Error: found %d left matches, but was expecting %d\n", matches, data.TrueMatchesLeft)
		}
		matches = CountRightMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
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
			matches = CountLeftMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
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
			matches = CountLeftMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
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
			matches = CountRightMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
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
			matches = CountRightMatches(data.SmallA, data.StartA, data.SmallB, data.StartB)
			if matches != data.TrueMatchesRight {
				b.Errorf("Error: wrong number of matches\n")
			}
		}
	}
}
