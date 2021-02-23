package interf

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"testing"
)

func TestBedLike(t *testing.T) {

	// bed reg  |----|
	// 1-base 1  2  3  4
	// 0-base 0  1  2  3
	// seq    A  C  G  T
	// If we have a vcf of 1 ACG -> A

	testBed := &bed.Bed{
		Chrom:      "chr1",
		ChromStart: 1,
		ChromEnd:   3}

	testVcf := &vcf.Vcf{
		Chr: "chr1",
		Pos: 1,
		Ref: "ACG",
		Alt: []string{"A"}}

	testAxt := &axt.Axt{
		RName:  "chr1",
		RStart: 2,
		REnd:   3}

	if testEquality(testBed, testVcf) && testEquality(testVcf, testAxt) {
		log.Println("BedLike is functioning properly")
	} else {
		t.Errorf("ERROR: Problem with BedLike methods")
	}

}

func testEquality(a BedLike, b BedLike) bool {
	testChr := a.GetChrom() == b.GetChrom()
	testStart := a.GetChromStart() == b.GetChromStart()
	testEnd := a.GetChromEnd() == b.GetChromEnd()

	if testChr && testStart && testEnd {
		return true
	} else {
		return false
	}
}
