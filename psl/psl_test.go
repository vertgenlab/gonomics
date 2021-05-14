package psl

import (
	"testing"
)

var psldata Psl = Psl{
	Match:       261,
	MisMatch:    54,
	RepeatMatch: 23,
	Ns:          0,
	QNumIns:     1,
	QBaseIns:    1,
	TNumIns:     0,
	TBaseIns:    0,
	Strand:      "+-",
	QName:       "NM_001001576",
	QSize:       376,
	QStart:      2,
	QEnd:        345,
	TName:       "chrI",
	TSize:       28185914,
	TStart:      8388245,
	TEnd:        8390701,
	BlockCount:  11,
	BlockSize:   []int{8, 17, 25, 29, 32, 29, 52, 36, 34, 44, 32},
	QList:       []int{2, 11, 28, 53, 84, 116, 145, 198, 235, 269, 313},
	TList:       []int{19795213, 19795237, 19795387, 19795565, 19795884, 19796147, 19796496, 19796787, 19797129, 19797357, 19797573},
}

func TestPslReader(t *testing.T) {
	readingFile := Read("testdata/pslLine.psl")
	if len(readingFile) != 1 {
		t.Errorf("Error: Please check file \"pslLine.psl\" and ensure it contains only 1 record...\n")
	} else {
		if !Equal(readingFile[0], psldata) {
			t.Errorf("Error: There was a parsing error that resulted in different psl structs...\n")
		}
	}
}

func BenchmarkPslReader(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Read("testdata/hgtBlast.psl")
	}
}
