package genePred

import (
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var genePred1 = GenePred{Id: "test", Chrom: "0", Strand: '+', TxStart: 0, TxEnd: 1001, CdsStart: 0, CdsEnd: 901, ExonNum: 4, ExonStarts: []int{0, 18, 500, 800}, ExonEnds: []int{3, 21, 503, 803}}
var genePred2 = GenePred{Id: "test", Chrom: "1", Strand: '+', TxStart: 0, TxEnd: 1001, CdsStart: 0, CdsEnd: 901, ExonNum: 4, ExonStarts: []int{0, 18, 500, 800}, ExonEnds: []int{3, 58, 602, 832}}
var genePreds = []GenePred{genePred1, genePred2}

var ReadTests = []struct {
	name string
	data []GenePred
}{
	{"testGenePred.gp", genePreds},
}

//func TestEqual(t *testing.T) {
//	for _, test := range ReadTests {
//		readingGenePred := Read(test.name)
//		checkEqual := Equal(genePreds[0], readingGenePred [0])
//		log.Print(checkEqual)
//	}
//}

func TestRead(t *testing.T) {
	for _, test := range ReadTests {
		actual := Read(test.name)
		equal, fields := AllAreEqual(test.data, actual)
		if !equal {
			t.Errorf("The %s file was not read correctly, or does not match genePred literals. Fields that did not match: %v", test.name, fields)
		}
	}
}

func TestWrite(t *testing.T) {
	for _, test := range ReadTests {
		Write("testWriting.gp", test.data)
	}
	fileio.EasyRemove("testWriting.gp")
}

//func TestCalcExonFrame(t *testing.T) {
//	for _, test := range ReadTests {
//		answer := CalcExonFrame(test.data[0])
//		log.Print(answer)
//	}
//}
