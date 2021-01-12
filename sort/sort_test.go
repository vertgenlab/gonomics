package sort

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"os"
	"testing"
)

func TestGirafExternalMergeSort(t *testing.T) {
	graph := simpleGraph.Read("testdata/mini.gg")
	sortOrder := simpleGraph.GetSortOrder(graph)
	GirafExternalMergeSort("testdata/miniReads.giraf", sortOrder, 100, "testSort.giraf")
	fileio.EasyRemove("testSort.giraf")
	fileio.EasyRemove("testSort.giraf.idx")
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
	fileio.EasyRemove("testdata/testWrite.idx")
}
