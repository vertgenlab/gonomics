package convert

import(
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

// ThreeWayFaToVcf takes in a three-way multiFa alignment and writes Vcf entries for segregating sites with the first
// entry as the reference and the last two fasta entries as the alt alleles.
// This is done by chromosome since a multiFa contains only one chromosome per file.
// This function only checks for substitutions, not indels.
func ThreeWayFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter) {
	var currRefPos, currAlnPos int = 0, 0

	if len(f) != 3 {
		log.Fatalf("ThreeWayFaToVcf expects a fasta input with three entries.")
	}

	for i := range f[0].Seq {
		if f[0].Seq[i] == dna.Gap || f[1].Seq[i] == dna.Gap || f[2].Seq[i] == dna.Gap {
			continue
		}
		if f[0].Seq[i] != f[1].Seq[i] || f[0].Seq[i] != f[2].Seq[i] { // normal substitution
			currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
			currAlnPos = i

			var altAllele []string
			var samples []vcf.Sample

			if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] == f[2].Seq[i] { // 1/0; substitution at allele1 only
				altAllele = []string{dna.BaseToString(f[1].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 0},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] == f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] { // 0/1; substitution at allele2 only
				altAllele = []string{dna.BaseToString(f[2].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{0, 1},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] && f[1].Seq[i] == f[2].Seq[i] { // 1/1; substitution at both allele1 and allele2, same alt
				altAllele = []string{dna.BaseToString(f[1].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 1},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			} else if f[0].Seq[i] != f[1].Seq[i] && f[0].Seq[i] != f[2].Seq[i] && f[1].Seq[i] != f[2].Seq[i] { // 1/2; substitution at both allele1 and allele2, different alt
				altAllele = []string{dna.BaseToString(f[1].Seq[i]), dna.BaseToString(f[2].Seq[i])}
				samples = []vcf.Sample{
					{
						Alleles:    []int16{1, 2},
						Phase:      []bool{false, false},
						FormatData: []string{""},
					},
				}
			}
			vcf.WriteVcf(out, vcf.Vcf{
				Chr:     chr,
				Pos:     currRefPos + 1, // 1-based indexing for VCF
				Id:      ".",
				Ref:     dna.BaseToString(f[0].Seq[i]), // ref allele in string format
				Alt:     altAllele,
				Qual:    100.0,
				Filter:  "PASS",
				Info:    ".",
				Format:  []string{"GT"},
				Samples: samples,
			})
		}
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

// PairwiseFaToVcf takes in a pairwise multiFa alignment and writes Vcf entries for segregating sites with the first
// entry as the reference and the second fasta entry as the alt allele.
// This will have to be done by chromosome, as a pairwise multiFa will only have two entries, thus containing one chromosome per file.
func PairwiseFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter, substitutionsOnly bool, retainN bool) {
	var pastStart, insertion, deletion bool = false, false, false //first bool checks to see if we have an insertion at the start of an alignment.
	var insertionAlnPos, deletionAlnPos int
	var currRefPos, currAlnPos int = 0, 0 //0 based, like fasta. Add 1 to get vcf pos.
	if len(f) != 2 {
		log.Fatalf("PairwiseFaToVcf expects a fasta input with two entries.")
	}

	for i := range f[0].Seq { //loop through alignment positions
		if f[0].Seq[i] == dna.Gap { //reference is gap (insertion)
			if pastStart {
				if !insertion {
					insertionAlnPos = i - 1
				}
				insertion = true
			}
		} else if f[0].Seq[i] != f[1].Seq[i] { // sequences diff at the same position
			pastStart = true
			if insertion { //catches the case where an insertion, now complete, is followed directly by a snp.
				if !substitutionsOnly {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], insertionAlnPos, currRefPos, currAlnPos)
					currAlnPos = insertionAlnPos //update currAlnPos
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
				}
			}
			if f[1].Seq[i] == dna.Gap { //alt is gap (deletion)
				if !deletion {
					deletionAlnPos = i - 1
				}
				deletion = true
			} else if deletion { //snp immediately follows the end of a deletion
				deletion = false
				if !substitutionsOnly {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], deletionAlnPos, currRefPos, currAlnPos)
					currAlnPos = deletionAlnPos
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
				}
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
						currAlnPos = i
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
					}
				} else {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
					currAlnPos = i
					vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //then add current diff
				}
			} else { //this case is for normal substitutions
				if f[0].Seq[i] == dna.N || f[1].Seq[i] == dna.N {
					if retainN {
						currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
						currAlnPos = i
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				} else {
					currRefPos = fasta.AlnPosToRefPosCounter(f[0], i, currRefPos, currAlnPos)
					currAlnPos = i
					if i < len(f[0].Seq)-1 { //if there is a next base to look at
						if f[0].Seq[i+1] != dna.Gap && f[1].Seq[i+1] != dna.Gap { //if neither alt nor ref is a gap in the next pos.
							vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
						} else if substitutionsOnly { //we would also write the subsitution in the case where we don't care if the substitution precedes an INDEL because we are only reporting substitutions
							vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
						}
						// Otherwise we won't report this substitution because it wil be part of the INDEL reporting.
					} else { //for a substitution in the final position, we need not check for subsequent INDELs
						vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[i]), Alt: []string{dna.BaseToString(f[1].Seq[i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
					}
				}
			}
			insertion = false
		} else if insertion { //case where ref and alt agree now but previous bases were part of an insertion.
			pastStart = true
			insertion = false
			if !substitutionsOnly {
				currRefPos = fasta.AlnPosToRefPosCounter(f[0], insertionAlnPos, currRefPos, currAlnPos)
				currAlnPos = insertionAlnPos
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BaseToString(f[0].Seq[insertionAlnPos]), Alt: []string{dna.BasesToString(f[1].Seq[insertionAlnPos:i])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}})
			}
		} else if deletion {
			pastStart = true
			deletion = false
			if !substitutionsOnly && deletionAlnPos >= 0 {
				//TODO: we do not currently save deletions if they occur at the start of an alignment.
				currRefPos = fasta.AlnPosToRefPosCounter(f[0], deletionAlnPos, currRefPos, currAlnPos)
				currAlnPos = deletionAlnPos
				vcf.WriteVcf(out, vcf.Vcf{Chr: chr, Pos: currRefPos + 1, Id: ".", Ref: dna.BasesToString(f[0].Seq[deletionAlnPos:i]), Alt: []string{dna.BaseToString(f[1].Seq[deletionAlnPos])}, Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}) //from deletion
			}
		}
	}
}

// NWayFaToVcf takes in a multiFa alignment with an arbitrary number of sequences and writes Vcf entries for segregating sites to out.
// This program treats the first entry as the reference and every subsequent sequence is treated as an independent haploid sequence, named after its fasta record name.
// Alignment columns containing a gap in any sequence are skipped. This function currently only supports substitutions.
func NWayFaToVcf(f []fasta.Fasta, chr string, out *fileio.EasyWriter) {
	var currRefPos, queryIdx, currAlnPos int = 0, 0, 0
	var refBase, queryBase dna.Base
	var isVariant, hasGap, found bool
	var alleleIdx int16

	if len(f) < 2 {
		log.Fatalf("convert.NWayFaToVcf expects a fasta input with at least two sequences.")
	}

	numQuery := len(f) - 1 // number of query sequences in the multiFa

	altIndex := make(map[dna.Base]int16, numQuery) //maps the alt base to its index (e.x. for a polymorphic Ref:A, Alt1:T, Alt2:C, map[C] = 2)
	altAlleleStrings := make([]string, 0, numQuery) // lists the observed alt alleles as strings, to be written out in the vcf. Reset each site
	samples := make([]vcf.Sample, numQuery) // to hold sample info for each vcf. Reset each site

	for alnPos := range f[0].Seq {
		isVariant = false
		hasGap = false
		clear(altIndex) // clearing the altIndex for each new site
		altAlleleStrings = altAlleleStrings[:0] // resetting the strings for each site

		for faIdx := range f {
			if f[faIdx].Seq[alnPos] == dna.Gap {
				hasGap = true
				break
			}
		}
		if hasGap {
			continue // skipping alignment columns with indels
		}

		refBase = f[0].Seq[alnPos]

		for queryIdx = 1; queryIdx <= numQuery; queryIdx++ {
			queryBase = f[queryIdx].Seq[alnPos]
			if queryBase != refBase {
				isVariant = true
				_, found = altIndex[queryBase]
				if !found {
					altAlleleStrings = append(altAlleleStrings, dna.BaseToString(queryBase))
					altIndex[queryBase] = int16(len(altAlleleStrings)) //this sets the altIndex 
				}
			}
		}

		if !isVariant { //skipping non-polymorphic sites
			continue
		}
		
		currRefPos = fasta.AlnPosToRefPosCounter(f[0], alnPos, currRefPos, currAlnPos)
		currAlnPos = alnPos

		for queryIdx = 1; queryIdx <= numQuery; queryIdx++ {
			queryBase = f[queryIdx].Seq[alnPos]
			alleleIdx = 0
			if queryBase != refBase {
				alleleIdx = altIndex[queryBase]
			}
			samples[queryIdx-1] = vcf.Sample {
				Alleles: []int16{alleleIdx},
				Phase: []bool{false},
				FormatData: []string{""},
			}
		}

		vcf.WriteVcf(out, vcf.Vcf{
			Chr:	chr,
			Pos:	currRefPos + 1, //1-based VCF
			Id:	".",
			Ref:	dna.BaseToString(refBase),
			Alt:	altAlleleStrings,
			Qual:	100.0,
			Filter:	"PASS",
			Info:	".",
			Format:	[]string{"GT"},
			Samples:	samples,
		})
	}
}
