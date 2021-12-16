package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"testing"
)

func TestGirafSort(t *testing.T) {
	graphFile := "../../sort/testdata/mini.gg"
	readsFile := "../../sort/testdata/miniReads.giraf"
	outFile := "miniReads_sorted_test.giraf"
	topoOrder := girafSort(readsFile, graphFile, 100, outFile)
	actual := giraf.Read(outFile)
	relabelIdBySortOrder(topoOrder, actual)

	var lastId uint32
	for i := range actual {
		if i == 0 {
			lastId = actual[i].Path.Nodes[0]
		}
		if actual[i].Path.Nodes[0] < lastId {
			t.Errorf("problem with girafSort")
		}
		lastId = actual[i].Path.Nodes[0]
	}

	fileio.EasyRemove(outFile)
	fileio.EasyRemove(outFile + ".idx")
}

func relabelIdBySortOrder(topoOrder []uint32, reads []*giraf.Giraf) {
	orderMap := make(map[uint32]uint32)
	for i := range topoOrder {
		orderMap[topoOrder[uint32(i)]] = uint32(i)
	}
	for i := range reads {
		for j := range reads[i].Path.Nodes {
			reads[i].Path.Nodes[j] = orderMap[reads[i].Path.Nodes[j]]
		}
	}
}
