package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/mHiC"
)

func writeInteractionSlice(toWrite []mHiC.Interaction, file *fileio.EasyWriter) {
	for i := range toWrite {
		fileio.WriteToFileHandle(file, toWrite[i].Full)
	}
}

func sameBins(a, b mHiC.Interaction) bool {
	if a.Chrom1 == b.Chrom1 && a.Chrom2 == b.Chrom2 && a.Bin1 == b.Bin1 && a.Bin2 == b.Bin2 {
		return true
	}
	return false
}

func filterInteractions(in, out string) {
	var currRead string
	var currInteractions []mHiC.Interaction
	var j, k int
	var found bool
	o := fileio.EasyCreate(out)
	interactChan := mHiC.GoReadInteractionToChan(in)
	for i := range interactChan {
		found = false
		if i.Uni {
			fileio.WriteToFileHandle(o, i.Full)
			continue
		}
		if i.ReadName != currRead {
			currRead = i.ReadName
			if k > 0 {
				writeInteractionSlice(currInteractions, o)
			}
			currInteractions = []mHiC.Interaction{i}
			k++
			continue
		}
		for j = range currInteractions {
			if sameBins(currInteractions[j], i) {
				found = true
			}
		}
		if !found {
			currInteractions = append(currInteractions, i)
		}
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func main() {
	filterInteractions("/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/TROPHOZOITES_PE.nameSort.validPairs", "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/TROPHOZOITES_PE.binFilter.validPairs")
}
