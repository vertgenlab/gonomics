package sort

import (
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"os"
	"testing"
)

func TestExternalMergeSort(t *testing.T) {
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
			log.Fatalln("ERROR: Problem with reading and writing sort indexes")
		}
	}

	os.Remove("testdata/testWrite")
}
