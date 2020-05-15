package axt

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
)

func HelperAxtToBed(a *Axt) *bed.Bed {
	return &bed.Bed{Chrom: a.RName, ChromStart: a.RStart - 1, ChromEnd: a.REnd - 1}
}

func NonOverlapToBed(axtfile []*Axt, bedfile []*bed.Bed, output string) {
	file := fileio.MustCreate(output)
	defer file.Close()
	for i := 0; i < len(bedfile); i++ {
		if NonOverlapBed(axtfile, bedfile[i]) {
			bed.WriteToFileHandle(file, bedfile[i], 5)
		}
	}
}

func NonOverlapBed(alignments []*Axt, peak *bed.Bed) bool {
	for _, a := range alignments {
		if OverlapAxtBed(a, peak) {
			return false
		}
	}
	return true
}
