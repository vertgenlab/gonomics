package chromInfo

import (
	"reflect"
	"testing"
)

func TestReading(t *testing.T) {
	correctChr1 := ChromInfo{"chr1", 249250621, 0}
	correctChr2 := ChromInfo{"chr2", 243199373, 1}
	correctChr3 := ChromInfo{"chr3", 198022430, 2}

	answerSlice := ReadToSlice("testdata/chromInfo.txt")
	if answerSlice[0] != correctChr1 ||
		answerSlice[1] != correctChr2 ||
		answerSlice[2] != correctChr3 {
		t.Errorf("ERROR: Problem reading chromInfo file to slice")
	}

	answerMap := ReadToMap("testdata/chromInfo.txt")

	if answerMap["chr1"] != correctChr1 ||
		answerMap["chr2"] != correctChr2 ||
		answerMap["chr3"] != correctChr3 {
		t.Errorf("ERROR: Problem reading chromInfo file to map")
	}

	if !reflect.DeepEqual(answerMap, SliceToMap(answerSlice)) {
		t.Errorf("ERROR: Problem converting chromInfo slice to map")
	}
}
