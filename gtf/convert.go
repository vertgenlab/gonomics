package gtf

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
)

// GeneToTssBed returns the positions of all TSSs from a Gene as a slice of single base-pair bed entries.
func GeneToTssBed(g Gene, c map[string]chromInfo.ChromInfo) []bed.Bed {
	return GeneToPromoterBed(g, c, 0 ,0)
}

//GenesToTssBed returns the position of all TSSs from a Gene map as a slice of single base-pair bed entries.
func GenesToTssBed(g map[string]*Gene, c map[string]chromInfo.ChromInfo) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	var currPromoters []bed.Bed
	var j int
	for _, i := range g {
		currPromoters = GeneToTssBed(*i, c)
		for j = range currPromoters {
			answer = append(answer, currPromoters[j])
		}
	}

	return answer
}

// GeneToPromoterBed produces a slice of beds from a gene containing the positions of promoters (TSS-500bp -> TSS+2kb)
// for all transcripts of the gene with the geneName in the Name field of the output Bed entries.
func GeneToPromoterBed(g Gene, c map[string]chromInfo.ChromInfo, upstream int, downstream int) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	var tBed bed.Bed
	for _, t := range g.Transcripts {
		if t.Strand {
			tBed = bed.Bed{Chrom: t.Chr, ChromStart:  numbers.Max(t.Start - upstream - 1, 0), ChromEnd: numbers.Min(t.Start + downstream, c[t.Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
		} else {
			tBed = bed.Bed{Chrom: t.Chr, ChromStart: numbers.Max(t.End - downstream - 1, 0), ChromEnd: numbers.Min(t.End + upstream, c[t.Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
		}
		answer = append(answer, tBed)
	}

	return answer
}

// GenesToPromoterBed produces a slice of beds from a set of genes containing the positions of all promoters
// for all transcripts for all genes with the geneID in the Name field of the output Bed entries.
func GenesToPromoterBed(g map[string]*Gene, c map[string]chromInfo.ChromInfo, upstream int, downstream int) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	var currPromoters []bed.Bed
	var j int

	for _, i := range g {
		currPromoters = GeneToPromoterBed(*i, c, upstream, downstream)

		for j = range currPromoters {
			answer = append(answer, currPromoters[j])
		}
	}

	return answer
}
