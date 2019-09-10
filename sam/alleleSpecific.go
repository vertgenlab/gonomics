package sam

import (
	//"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"

	"strings"
	//"fmt"
)

//Need to filter VCF for only single base SNPS
func AlleleExpression(v []*vcf.Vcf, chromSize map[string]*chromInfo.ChromInfo) (map[int64][]dna.Base, map[int64][]dna.Base) {
	ref := make(map[int64][]dna.Base)
	alt := make(map[int64][]dna.Base)
	var curr *vcf.Vcf
	var faIdx int64
	for i := 0; i < len(v); i++ {
		curr = v[i]
		faIdx = chromSize[curr.Chr].Order << 32

		//convert to zero based as standard for programing
		ref[faIdx|curr.Pos-1] = dna.StringToBases(curr.Ref)
		alt[faIdx|curr.Pos-1] = dna.StringToBases(curr.Alt)
	}
	return ref, alt
}

func SnpSort(align []*SamAln, v []*vcf.Vcf, chromSize map[string]*chromInfo.ChromInfo) ([]*SamAln, []*SamAln) {
	var fresh []*SamAln
	var marine []*SamAln
	var fMap, mMap = AlleleExpression(v, chromSize)
	//fmt.Println("length of fMap", len(v))
	//fmt.Println("length of SAM", len(align))
	var i, j, fSnpCount, mSnpCount int
	var currRefPos, currQueryPos, k, faIdx int64

	var numDiscardReads int64
	for i = 0; i < len(align); i++ {
		fSnpCount, mSnpCount = 0, 0
		//convert sam (1-based) to zero base
		currRefPos = align[i].Pos - 1
		currQueryPos = 0
		faIdx = chromSize[align[i].RName].Order << 32

		for j = 0; j < len(align[i].Cigar); j++ {
			switch align[i].Cigar[j].Op {
			case 'S':
				currRefPos += align[i].Cigar[j].RunLength - 1
				currQueryPos += align[i].Cigar[j].RunLength - 1
			case 'I':
				//ask craig
				currQueryPos += align[i].Cigar[j].RunLength - 1
			case 'D':
				//ask craig
				currRefPos += align[i].Cigar[j].RunLength - 1
			case 'M':
				//check for mismatches
				for k = 0; k < align[i].Cigar[j].RunLength; k++ {
					//fmt.Println("Ref: ", dna.BasesToString(fMap[faIdx|currRefPos+k]), "Alt: ", dna.BasesToString(mMap[faIdx|currRefPos+k]))
					if strings.Compare(string(align[i].Seq[currQueryPos+k]), dna.BasesToString(fMap[faIdx|currRefPos+k])) == 0 {
						fSnpCount++
					} else if strings.Compare(string(align[i].Seq[currQueryPos+k]), dna.BasesToString(mMap[faIdx|currRefPos+k])) == 0 {
						mSnpCount++
					} else {
						//do nothing, this includes if key (position) does not exist
					}
				}
				currRefPos += align[i].Cigar[j].RunLength - 1
				currQueryPos += align[i].Cigar[j].RunLength - 1
			}
		}
		if fSnpCount > mSnpCount {
			fresh = append(fresh, align[i])
		} else if mSnpCount > fSnpCount {
			marine = append(marine, align[i])
		} else {
			numDiscardReads++
		}
	}
	return fresh, marine
}
