package gtf

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
)

// GeneToTssBed returns the positions of all TSSs from a Gene as a slice of single base-pair bed entries.
func GeneToTssBed(g Gene, c map[string]chromInfo.ChromInfo) []bed.Bed {
	return GeneToPromoterBed(g, c, 0, 0)
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

// GenesToCanonicalTrancriptsTssBed turns an input map of [geneId]*Gene structs, finds the canonical
// transcript (defined as the longest coding sequence), and turns the TSS of this trancript into a Bed
//struct.
func GenesToCanonicalTranscriptsTssBed(g map[string]*Gene, c map[string]chromInfo.ChromInfo) []bed.Bed {
	var answer []bed.Bed
	for _, i := range g {
		answer = append(answer, GeneToCanonicalTssBed(*i, c))
	}
	return answer
}

// GeneToCanonicalTssBed converts a single Gene struct into a Bed representing the TSS position of
// the canonical transcript.
func GeneToCanonicalTssBed(g Gene, c map[string]chromInfo.ChromInfo) bed.Bed {
	return GeneToCanonicalBed(g, c, 0, 0)
}

// GenesToCanonicalBeds converts all genes in a map[string]*Gene to a []bed.Bed, where each bed
// represeents the promoter region of the canonical transcript, defined by user-specified upstream
// and downstream distances from the TSS.
func GenesToCanonicalBeds(g map[string]*Gene, c map[string]chromInfo.ChromInfo, upstream int, downstream int) []bed.Bed {
	var answer []bed.Bed
	for _, i := range g {
		answer = append(answer, GeneToCanonicalBed(*i, c, upstream, downstream))
	}
	return answer
}

// GeneToCanonicalBed converts a Gene struct into a bed representing the promoter region of the
// canonical transcript. The user species the bases upstream and downstream of the TSS
// which will define the promoter region.
func GeneToCanonicalBed(g Gene, c map[string]chromInfo.ChromInfo, upstream int, downstream int) bed.Bed {
	MoveCanonicalToZero(&g)
	if g.Transcripts[0].Strand {
		return bed.Bed{Chrom: g.Transcripts[0].Chr, ChromStart: numbers.Max(g.Transcripts[0].Start-upstream-1, 0), ChromEnd: numbers.Min(g.Transcripts[0].Start+downstream, c[g.Transcripts[0].Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
	} else {
		return bed.Bed{Chrom: g.Transcripts[0].Chr, ChromStart: numbers.Max(g.Transcripts[0].End-downstream-1, 0), ChromEnd: numbers.Min(g.Transcripts[0].End+upstream, c[g.Transcripts[0].Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
	}
}

// GeneToPromoterBed produces a slice of beds from a gene containing the positions of promoters (TSS-500bp -> TSS+2kb)
// for all transcripts of the gene with the geneName in the Name field of the output Bed entries.
func GeneToPromoterBed(g Gene, c map[string]chromInfo.ChromInfo, upstream int, downstream int) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	var tBed bed.Bed
	for _, t := range g.Transcripts {
		if t.Strand {
			tBed = bed.Bed{Chrom: t.Chr, ChromStart: numbers.Max(t.Start-upstream-1, 0), ChromEnd: numbers.Min(t.Start+downstream, c[t.Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
		} else {
			tBed = bed.Bed{Chrom: t.Chr, ChromStart: numbers.Max(t.End-downstream-1, 0), ChromEnd: numbers.Min(t.End+upstream, c[t.Chr].Size), Name: g.GeneName, FieldsInitialized: 4}
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
