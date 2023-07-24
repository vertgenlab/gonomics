package gtf

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/vcf"
)

// VariantArrayOverlap takes a vcf record and a map of bool slices (chrom name maps to a bool for each base in that chrom).
// The bool slice encodes the presense/absense of some genomic feature and true is returned if the vcf position overlaps that feature.
func VariantArrayOverlap(v *vcf.Vcf, a map[string][]bool) bool {
	return a[v.Chr][v.Pos-1]
}

// ExonBoolArray returns a map of chromosome names to bool slices.  The bool is true if that position lies in an exon.
// A map of chromosome name to the information for that chromosome is needed to know the length of the retuned bool slices.
func ExonBoolArray(g map[string]*Gene, c map[string]*chromInfo.ChromInfo) map[string][]bool {
	var answer map[string][]bool
	for k := range c {
		answer[k] = make([]bool, c[k].Size)
	}

	for k := range g {
		for i := 0; i < len(g[k].Transcripts); i++ {
			for j := 0; j < len(g[k].Transcripts[i].Exons); j++ {
				//End is not minus 1 because the 1 basing is canceled out by the closed right interval of gtf.
				for m := g[k].Transcripts[i].Exons[j].Start - 1; m < g[k].Transcripts[i].Exons[j].End; m++ {
					answer[k][m] = true
				}
			}
		}
	}
	return answer
}

// CdsBoolArray returns a map of chromosome names to bool slices.  The bool is true if that position lies in a cds (protein-coding) region.
// A map of chromosome name to the information for that chromosome is needed to know the length of the retuned bool slices.
func CdsBoolArray(g map[string]*Gene, c map[string]*chromInfo.ChromInfo) map[string][]bool {
	var answer map[string][]bool
	for k := range c {
		answer[k] = make([]bool, c[k].Size)
	}

	for k := range g {
		for i := 0; i < len(g[k].Transcripts); i++ {
			for j := 0; j < len(g[k].Transcripts[i].Exons); j++ {
				if g[k].Transcripts[i].Exons[j].Cds != nil {
					for m := g[k].Transcripts[i].Exons[j].Cds.Start - 1; m < g[k].Transcripts[i].Exons[j].Cds.End; m++ {
						answer[k][m] = true
					}
				}
			}
		}
	}
	return answer
}

// FiveUtrBoolArray returns a map of chromosome names to bool slices.  The bool is true if that position lies in a 5' UTR.
// A map of chromosome name to the information for that chromosome is needed to know the length of the retuned bool slices.
func FiveUtrBoolArray(g map[string]*Gene, c map[string]*chromInfo.ChromInfo) map[string][]bool {
	var answer map[string][]bool
	for k := range c {
		answer[k] = make([]bool, c[k].Size)
	}

	for k := range g {
		for i := 0; i < len(g[k].Transcripts); i++ {
			for j := 0; j < len(g[k].Transcripts[i].Exons); j++ {
				if g[k].Transcripts[i].Exons[j].FiveUtr != nil {
					for m := g[k].Transcripts[i].Exons[j].FiveUtr.Start - 1; m < g[k].Transcripts[i].Exons[j].FiveUtr.End; m++ {
						answer[k][m] = true
					}
				}
			}
		}
	}
	return answer
}

// ThreeUtrBoolArray returns a map of chromosome names to bool slices.  The bool is true if that position lies in a 3' UTR.
// A map of chromosome name to the information for that chromosome is needed to know the length of the retuned bool slices.
func ThreeUtrBoolArray(g map[string]*Gene, c map[string]*chromInfo.ChromInfo) map[string][]bool {
	var answer map[string][]bool
	for k := range c {
		answer[k] = make([]bool, c[k].Size)
	}

	for k := range g {
		for i := 0; i < len(g[k].Transcripts); i++ {
			for j := 0; j < len(g[k].Transcripts[i].Exons); j++ {
				if g[k].Transcripts[i].Exons[j].ThreeUtr != nil {
					for m := g[k].Transcripts[i].Exons[j].ThreeUtr.Start - 1; m < g[k].Transcripts[i].Exons[j].ThreeUtr.End; m++ {
						answer[k][m] = true
					}
				}
			}
		}
	}
	return answer
}
