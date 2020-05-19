package bed

import (
	"github.com/vertgenlab/gonomics/axt"
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

	testBed := &Bed{
		Chrom: "chr1",
		ChromStart: 1,
		ChromEnd: 3}

	testVcf := &vcf.Vcf{
		Chr: "chr1",
		Pos: 1,
		Ref: "ACG",
		Alt: "A"}

	testAxt := &axt.Axt{
		RName: "chr1",
		RStart: 2,
		REnd: 3}

	if testEquality(testBed, testVcf) && testEquality(testVcf, testAxt) {
		log.Println("BedLike is functioning properly")
	} else {
		t.Errorf("ERROR: Problem with BedLike methods")
	}

}

func testEquality(a BedLike, b BedLike) bool {
	testChr := a.GetChr() == b.GetChr()
	testStart := a.GetStart() == b.GetStart()
	testEnd := a.GetEnd() == b.GetEnd()

	if testChr && testStart && testEnd {
		return true
	} else {
		return false
	}
}
