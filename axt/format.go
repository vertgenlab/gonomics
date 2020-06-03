package axt

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
)

func AxtBed(aRec *Axt, checkTarget bool) *bed.Bed {
	if checkTarget {
		return &bed.Bed{Chrom: aRec.RName, ChromStart: aRec.RStart - 1, ChromEnd: aRec.REnd - 1, Name: aRec.QName, Score: int64(aRec.Score)}
	} else {
		return queryAxtBed(aRec)
	}
}

//I believe target is positive, and query and be positive strand or negatiove strand
func queryAxtBed(axtRec *Axt) *bed.Bed {
	if axtRec.QStrandPos {
		return &bed.Bed{Chrom: axtRec.QName, ChromStart: axtRec.QStart - 1, ChromEnd: axtRec.QEnd - 1, Name: axtRec.QName, Score: int64(axtRec.Score)}
	} else {
		return &bed.Bed{Chrom: axtRec.QName, ChromStart: int64(len(axtRec.QSeq)) - axtRec.QEnd - 1, ChromEnd: int64(len(axtRec.QSeq)) - axtRec.QStart - 1, Name: axtRec.QName, Score: int64(axtRec.Score)}
	}
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
