package gtf

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/vcf"
)

func VariantArrayOverlap(v *vcf.Vcf, a map[string][]bool) bool {
	return a[v.Chr][v.Pos-1]
}

// returns a map of chromosome names to bool arrays, true if that position lies in an exon.
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
