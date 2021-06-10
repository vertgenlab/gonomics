package bed

import (
	"log"
	"testing"
)

var r1 Bed = Bed{Chrom: "chr1", ChromStart: 1, ChromEnd: 2}
var r2 Bed = Bed{Chrom: "chr1", ChromStart: 4, ChromEnd: 5}
var r3 Bed = Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 11}
var r4 Bed = Bed{Chrom: "chr1", ChromStart: 12, ChromEnd: 14}

var querySlice []*Bed = []*Bed{&r1, &r2, &r3, &r4}
var chromSize int = 15

var i1 Bed = Bed{Chrom: "chr1", ChromStart: 0, ChromEnd: 1}
var i2 Bed = Bed{Chrom: "chr1", ChromStart: 2, ChromEnd: 4}
var i3 Bed = Bed{Chrom: "chr1", ChromStart: 5, ChromEnd: 10}
var i4 Bed = Bed{Chrom: "chr1", ChromStart: 11, ChromEnd: 12}
var i5 Bed = Bed{Chrom: "chr1", ChromStart: 14, ChromEnd: 15}
var ans []*Bed = []*Bed{&i1, &i2, &i3, &i4, &i5}

func TestInvertBedRegions(t *testing.T) {
	testSlice := InvertRegions(querySlice, chromSize)
	//log.Printf("len=%d\n", len(testSl))
	for i, b := range testSlice {
		if !Equal(b, ans[i]) {
			log.Printf("Bed=%s\n", ToString(b, 3))
		}
	}
}
