package axt

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
)

func FilterOverlap(axtfile []*Axt, bedfile []*bed.Bed, output string, target bool, query bool) {
	var i, j int
	genome := make(map[int]*Axt)
	var key bool
	for i = 0; i < len(axtfile); i++ {
		if !ContainsNs(axtfile[i], target, query) {
			for j = 0; j < len(bedfile); j++ {
				if OverlapAxtBed(axtfile[i], bedfile[j]) {
					_, key = genome[i]
					if !key {
						genome[i] = axtfile[i]
					}
				}
			}
		}
	}
	file := fileio.EasyCreate(output)
	defer file.Close()
	for i = 0; i < len(genome); i++ {
		WriteToFileHandle(file, genome[i], i)
	}
}

func NonOverlapFilter(axtfile []*Axt, bedfile []*bed.Bed, output string, target bool, query bool) {
	file := fileio.EasyCreate(output)
	defer file.Close()
	var index int = 0
	var i int
	for i = 0; i < len(axtfile); i++ {
		if !ContainsNs(axtfile[i], target, query) {
			if NonOverlapAxtBed(axtfile[i], bedfile) {
				WriteToFileHandle(file, axtfile[i], index)
				index++
			}
		}
	}
}

func NonOverlapAxtBed(alignment *Axt, peaks []*bed.Bed) bool {
	for _, b := range peaks {
		if OverlapAxtBed(alignment, b) {
			return false
		}
	}
	return true
}

func OverlapAxtBed(alpha *Axt, beta *bed.Bed) bool {
	if (common.MaxInt64(alpha.RStart-1, beta.ChromStart) < common.MinInt64(alpha.REnd-1, beta.ChromEnd)) && strings.Compare(alpha.RName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func AxtVcfOverlap(alpha *Axt, beta *vcf.Vcf) bool {
	if common.MaxInt64(alpha.RStart-1, beta.Pos-1) < common.MinInt64(alpha.REnd-1, int64(len(dna.StringToBases(beta.Ref)))) && strings.Compare(alpha.RName, beta.Chr) == 0 {
		return true
	} else {
		return false
	}
}

func ContainsNs(block *Axt, target bool, query bool) bool {
	if target && query {
		if dna.CountBaseInterval(block.RSeq, dna.N, 0, len(block.RSeq)) != 0 || dna.CountBaseInterval(block.QSeq, dna.N,
			0, len(block.QSeq)) != 0 {
			return true
		}
	}
	if target && !query {
		if dna.CountBaseInterval(block.RSeq, dna.N, 0, len(block.RSeq)) != 0 {
			return true
		}
	}
	if query && !target {
		if dna.CountBaseInterval(block.QSeq, dna.N, 0, len(block.QSeq)) != 0 {
			return true
		}
	}
	return false
}
