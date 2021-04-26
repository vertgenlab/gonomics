package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"strings"
)

//AppendAncestor adds the ancestral allele state (defined by input bases) to the INFO column of a vcf entry.
func AppendAncestor(g Vcf, b []dna.Base) Vcf {
	if g.Info == "." {
		g.Info = fmt.Sprintf("AA=%s", dna.BasesToString(b))
	} else {
		g.Info = fmt.Sprintf("%s;AA=%s", g.Info, dna.BasesToString(b))
	}
	return g
}

//IsRefAncestor returns true if the reference allele in the record matches the ancestral allele in the Info annotation, false otherwise.
func IsRefAncestor(g Vcf) bool {
	return dna.BasesToString(QueryAncestor(g)) == g.Ref
}

//IsAltAncestor returns true if the first alt allele in the record matches the ancestral allele in the Info annotation, false otherwise.
func IsAltAncestor(g Vcf) bool {
	return dna.BasesToString(QueryAncestor(g)) == g.Alt[0]
}

//QueryAncestor finds the AA INFO from a VCF struct and returns the base of the ancestral allele.
func QueryAncestor(g Vcf) []dna.Base {
	if g.Info == "." {
		return nil //or should this log.Fatalf out? Depends on whether we have vcf with partial annotation
	}
	fields := strings.Split(g.Info, ";")
	var f []string
	for i := 0; i < len(fields); i++ {
		f = strings.Split(fields[i], "=")
		if f[0] == "AA" {
			return dna.StringToBases(f[1])
		}
	}
	return nil
}

//HasAncestor returns true if a VCF record is annotated with an ancestor allele in the Info column, false otherwise.
func HasAncestor(g Vcf) bool {
	return QueryAncestor(g) != nil
}

//AnnotateAncestorFromMultiFa adds the ancestral state to a VCF variant by inspecting a pairwise fasta of the reference genome and an ancestor sequence.
//records is a pairwise multiFa where the first entry is the reference genome and the second entry is the ancestor.
func AnnotateAncestorFromMultiFa(g Vcf, records []fasta.Fasta) Vcf {
	p := fasta.RefPosToAlnPos(records[0], int(g.Pos)-1) //get the alignment position of the variant
	//DEBUG: fmt.Printf("RefSeq: %s\n", dna.BasesToString(records[0].Seq))
	//DEBUG: fmt.Printf("Alignment pos: %v. Base at p: %s. Base at p+1: %s.\n", p, dna.BaseToString(records[0].Seq[p]), dna.BaseToString(records[0].Seq[p+1]))
	var AncestralAllele []dna.Base
	var insertionEnd int
	if records[0].Seq[p+1] == dna.Gap { //true in the case of insertions, as there is a gap in the reference after the variant position.
		insertionEnd = p + 1
		//DEBUG: fmt.Printf("Found insertion.\n")
		for insertionEnd < len(records[0].Seq) {
			if records[0].Seq[insertionEnd] != dna.Gap {
				break
			}
			insertionEnd++
		}
		AncestralAllele = records[1].Seq[p:insertionEnd]
	} else { //No gaps after the pos in ref means we have a deletion or snp, so the ancestral allele is then just the base at p.
		AncestralAllele = records[1].Seq[p : p+1]
	}
	g = AppendAncestor(g, AncestralAllele)
	return g
}

//AncestorFlagToHeader adds an ##INFO line to a vcfHeader to include information about the AA flag for ancestral alleles.
func AncestorFlagToHeader(h Header) Header {
	var lastInfoIndex int
	var firstFormatIndex int = -1
	var seenInfo, seenFormat, firstTime bool
	var AncestorLine string = "##INFO=<ID=AA,Number=1,Type=String,Description=\"AncestralAllele\">"
	var AncestorLineSlice []string = make([]string, 1)
	AncestorLineSlice[0] = AncestorLine
	firstTime = true

	for i := 0; i < len(h.Text); i++ {
		fields := strings.Split(h.Text[i], "=")
		if fields[0] == "##INFO" {
			if !seenInfo {
				seenInfo = true
			}
		} else if seenInfo && firstTime {
			lastInfoIndex = i
			firstTime = false
		}
		if fields[0] == "##FORMAT" && !seenFormat {
			seenFormat = true
			firstFormatIndex = i
		}
	}
	//if we didn't see info columns, we append this above the first FORMAT line
	if !seenInfo {
		if firstFormatIndex == -1 { //in this case we didn't see format lines or info lines, so the AA flag is simply appended to the header
			h.Text = append(h.Text, AncestorLine)
		}
		h.Text = append(h.Text[:firstFormatIndex], append(AncestorLineSlice, h.Text[firstFormatIndex:]...)...)
	} else { //otherwise, we insert the new headerline after the last info line.
		//DEBUG: fmt.Printf("Length of header: %v. LastInfoIndex: %v.\n", len(h.Text), lastInfoIndex)
		h.Text = append(h.Text[:lastInfoIndex], append(AncestorLineSlice, h.Text[lastInfoIndex:]...)...)
	}
	return h
}
