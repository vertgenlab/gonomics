package sort

import (
	"fmt"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"os"
	"testing"
	"unsafe"
)

type ByteTest struct {
	d	int32
	s 	int32
	a 	byte
	c 	byte
}

func TestExternalMergeSort(t *testing.T) {
	var test ByteTest
	val := unsafe.Sizeof(test)
	fmt.Println(val)
	graph := simpleGraph.Read("testdata/mini.gg")
	sortOrder := simpleGraph.GetSortOrder(graph)
	ExternalMergeSort("testdata/miniReads.giraf", sortOrder, 100, "testSort.giraf")
}

func TestReadAndWriteIdx(t *testing.T) {
	inputOrder := []uint32{0, 5, 3, 2, 0}
	writeIdx("testdata/testWrite", inputOrder)
	ouputOrder := ReadIdx("testdata/testWrite.idx")

	for i := 0; i < len(inputOrder); i++ {
		if inputOrder[i] != ouputOrder[i] {
			t.Error("ERROR: Problem with reading and writing sort indexes")
		}
	}

	os.Remove("testdata/testWrite")
}
