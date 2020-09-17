package binaryGiraf

import (
	"github.com/vertgenlab/gonomics/giraf"
	"reflect"
	"sync"
	"testing"
)

const (
	TestGiraf   = "testdata/test.giraf"
	TestGraph   = "testdata/test.gg"
	TestGiraffe = "testdata/test.giraf.fe"
)

func TestReadToSlice(t *testing.T) {
	data := Read(TestGiraffe, TestGraph)
	correct := giraf.Read(TestGiraf)
	if len(data) != len(correct) || len(data) != 1 {
		t.Errorf("func Read: Decompressed giraf files are not correct length")
	}
	if !reflect.DeepEqual(data[0], *(correct[0])) {
		t.Errorf("func Read: Decompressed giraf file is not equal to uncompressed giraf file")
	}
}

func TestReadToChan(t *testing.T) {
	data := GoReadToChan(TestGiraffe, TestGraph)
	correct := giraf.Read(TestGiraf)

	for record := range data {
		if !reflect.DeepEqual(record, *(correct[0])) {
			t.Errorf("func ReadToChan: Decompressed giraf file is not equal to uncompressed giraf file")
		}
	}
}

func TestGirafChanToBinary(t *testing.T) {
	data := giraf.GoReadToChan(TestGiraf)
	correct := giraf.Read(TestGiraf)
	var wg sync.WaitGroup
	wg.Add(1)
	go GirafChanToBinary("testdata/tmp.giraf.fe", data, &wg)
	wg.Wait()

	answer := Read("testdata/tmp.giraf.fe", TestGraph)
	if len(answer) != len(correct) || len(answer) != 1 {
		t.Errorf("func GirafChanToBinary: Decompressed giraf files are not correct length")
	}
	if !reflect.DeepEqual(answer[0], *(correct[0])) {
		t.Errorf("func GirafChanToBinary: Written giraf.fe does not equal uncompressed giraf")
	}
}
