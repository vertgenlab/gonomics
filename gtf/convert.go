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

// GenesToTssBed returns the position of all TSSs from a Gene map as a slice of single base-pair bed entries.
func GenesToTssBed(g map[string]*Gene, c map[string]chromInfo.ChromInfo, merge bool) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	var currPromoters []bed.Bed
	for _, i := range g {
		currPromoters = GeneToTssBed(*i, c)
		answer = append(answer, currPromoters...)
	}

	if merge {
		answer = bed.MergeBeds(answer)
	}

	return answer
}

// GenesToCanonicalTranscriptsTssBed turns an input map of [geneId]*Gene structs, finds the canonical
// transcript (defined as the longest coding sequence), and turns the TSS of this trancript into a Bed
// struct.
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
// represents the promoter region of the canonical transcript, defined by user-specified upstream
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

	for _, i := range g {
		currPromoters = GeneToPromoterBed(*i, c, upstream, downstream)
		answer = append(answer, currPromoters...)
	}

	return answer
}

// GenesToBedFirstTwoCodonBases takes a map[string]*Gene and returns a []bed.Bed containing the first two bases
// of each codon. The result will be returned in coordinate sorted order.
func GenesToBedFirstTwoCodonBases(genes map[string]*Gene) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	for currGene := range genes {
		answer = append(answer, geneToBedFirstTwoCodonBases(*genes[currGene])...)
	}
	bed.SortByCoord(answer)
	return answer
}

// geneToBedFirstTwoCodonBases is a helper function of GenesToBedFirstTwoCodonBases. It takes a single Gene
// and produces a slice of bed.Bed structs representing the first two bases of each codon for every coding exon
// of every transcript.
func geneToBedFirstTwoCodonBases(g Gene) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	for _, t := range g.Transcripts {
		for _, e := range t.Exons {
			if e.Cds != nil {
				answer = append(answer, cdsToBedsFirstTwoCodonBases(*e.Cds, t.Chr, t.Strand)...)
			}
		}
	}
	return answer
}

// cdsToBedsFirstTwoCodonBases is a helper function of geneToBedFirstTwoCodonBases. It takes a Cds struct and
// returns a slice of bed.Bed structs representing the first two bases of each codon. This function is aware of
// both frame and strand. Useful for conservation analysis as this excludes degenerate sites.
func cdsToBedsFirstTwoCodonBases(c Cds, chrom string, strand bool) []bed.Bed {
	var answer = make([]bed.Bed, 0)
	// we loop through all complete codons in the cds (ranging from start + frame - 1).
	// the -1 handles the fact that GTF is one-based and BED is 0-based. The + 2 for the end
	// reflects the fact that BED includes the start and excludes the end.
	for currPos := c.Start + c.Frame - 1; currPos+2 <= c.End; currPos += 3 {
		if strand {
			answer = append(answer, bed.Bed{Chrom: chrom, ChromStart: currPos, ChromEnd: currPos + 2, FieldsInitialized: 3})
		} else {
			answer = append(answer, bed.Bed{Chrom: chrom, ChromStart: currPos + 1, ChromEnd: currPos + 3, FieldsInitialized: 3})
		}
	}
	return answer
}
