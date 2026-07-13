// Notes:
// set divBedDir = "/work/yl726/pfaDivBedOutputs", files like chr5.107183644.107184616.Consensus_HAQER_0316.bed
// set vcfDir    = "/hpc/group/vertgenlab/raven/PrimateT2T_20way/selection_hs1/populationData/Filtered", files like chr22.segSite.vcf.gz
// (only for reference now, not used in program) divSegDir = /work/yl726/divVcfMatch, files like Consensus_HAQER_0006.divSeg.tsv, example divSeg line: chr10	12130633	chr10_12130633_G_A	G	A	255	.	.	GT	0|0	chr10	12130632	12130633	Consensus_HAQER_0006	0	.	G	A
// Run this beforehand: cat /work/yl726/pfaDivBedOutputs/*.bed > allPfaDivs.bed
// CLI below
// divVcfMatch -div allPfaDivs.bed -vcfDir /hpc/group/vertgenlab/raven/PrimateT2T_20way/selection_hs1/populationData/Filtered -out allDivSegs.tsv

package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/vcf"
)

// DivBed represents one divergence line read from the BED-style input.
// We expect BED lines like:
// chrom start end name score strand div anc
type DivBed struct {
	Chrom   string
	Start0  int // 0-based start from BED
	End0    int
	Name    string // HAQER id (col 4)
	DivBase string // col 7
	AncBase string // col 8
}

// DivMap: chrom -> pos1 -> []DivBed
type DivMap map[string]map[int][]DivBed

// readDivBed reads the input file produced by your pfaDivBed outputs (concatenated).
// It expects each line to be a BED-like line with at least 8 tab-separated columns:
// chrom (1), start (2), end (3), name (4), score (5), strand (6), div (7), anc (8)
func readDivBed(path string) (DivMap, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	ans := make(DivMap)
	scanner := bufio.NewScanner(f)
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := scanner.Text()
		if strings.TrimSpace(line) == "" {
			continue
		}
		// Accept either tab-separated or whitespace
		fields := strings.Fields(line)
		// But try splitting by tab first to preserve empty fields (if any)
		if strings.Contains(line, "\t") {
			fields = strings.Split(line, "\t")
		}
		if len(fields) < 8 {
			return nil, fmt.Errorf("line %d: expected >=8 columns (chrom start end name score strand div anc), got %d: %q", lineNo, len(fields), line)
		}
		chrom := fields[0]
		start0, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid start: %v", lineNo, err)
		}
		end0, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid end: %v", lineNo, err)
		}
		name := fields[3]
		// fields[4] score, fields[5] strand are ignored by the program but kept in input
		div := fields[6]
		anc := strings.ToUpper(fields[7])

		rec := DivBed{
			Chrom:   chrom,
			Start0:  start0,
			End0:    end0,
			Name:    name,
			DivBase: div,
			AncBase: anc,
		}

		pos1 := start0 + 1 // VCF 1-based pos

		if _, ok := ans[chrom]; !ok {
			ans[chrom] = make(map[int][]DivBed)
		}
		ans[chrom][pos1] = append(ans[chrom][pos1], rec)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return ans, nil
}

// alleleCountsForVcf computes AC and AN for desired allele index in vcf.Vcf
// desiredAlleleIdx: 0 => REF, 1.. => ALT index+1
func alleleCountsForVcf(v vcf.Vcf, desiredAlleleIdx int) (ac int, an int) {
	ac = 0
	an = 0

	for _, s := range v.Samples {
		// prefer pre-parsed Alleles if present
		if len(s.Alleles) > 0 {
			for _, a := range s.Alleles {
				if int(a) >= 0 { // allele index -1 indicates missing
					an++
					if int(a) == desiredAlleleIdx {
						ac++
					}
				}
			}
			continue
		}
		// fallback: parse GT-like field from FormatData[0]
		if len(s.FormatData) > 0 {
			gt := s.FormatData[0]
			if gt == "." || gt == "./." || gt == ".|." {
				continue
			}
			alleles := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
			for _, aStr := range alleles {
				if aStr == "." {
					continue
				}
				idx, err := strconv.Atoi(aStr)
				if err != nil {
					continue
				}
				an++
				if idx == desiredAlleleIdx {
					ac++
				}
			}
			continue
		}
		// no genotype info: skip
	}
	return ac, an
}

func processUsingVcfDir(divs DivMap, vcfDir string, minFreq float64, requireSnp bool, out *os.File) {
	for chrom, chromMap := range divs {
		vcfName := chrom + ".segSite.vcf.gz"
		vcfPath := filepath.Join(vcfDir, vcfName)

		if _, err := os.Stat(vcfPath); os.IsNotExist(err) {
			log.Printf("warning: VCF for %s not found at %s — skipping\n", chrom, vcfPath)
			continue
		}

		vcfChan, _ := vcf.GoReadToChan(vcfPath)
		for rec := range vcfChan {
			pos := rec.Pos
			divList, ok := chromMap[pos]
			if !ok {
				continue
			}
			// handle each divergence record at this pos
			for _, d := range divList {
				anc := d.AncBase
				if requireSnp && len(rec.Ref) != 1 {
					continue
				}
				desiredAlleleIdx := -1
				matchedAlt := "." // if anc matched an ALT, which one?
				if anc == rec.Ref {
					desiredAlleleIdx = 0
				} else {
					for ai, alt := range rec.Alt {
						if strings.ToUpper(alt) == anc {
							desiredAlleleIdx = ai + 1
							matchedAlt = alt
							break
						}
					}
				}
				if desiredAlleleIdx == -1 {
					// ancestor base not present as REF/ALT in this VCF record; skip
					continue
				}

				ac, an := alleleCountsForVcf(rec, desiredAlleleIdx)
				if an == 0 {
					continue
				}
				af := float64(ac) / float64(an)
				if af > minFreq {
					// Output: chr start(bed) vcfId vcfRefBase vcfAltBase vcfFrequencyOfHaqerAncBase haqerId haqerDivBase haqerAncBase
					vcfAltOut := matchedAlt
					if anc == rec.Ref {
						vcfAltOut = "." // anc is ref; no specific matching ALT
					}
					start0 := d.Start0
					fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\n",
						d.Chrom, start0, rec.Id, rec.Ref, vcfAltOut, af, d.Name, d.DivBase, d.AncBase)
				}
			}
		}
	}
}

func processUsingSingleVcf(divs DivMap, vcfPath string, minFreq float64, requireSnp bool, out *os.File) {
	vcfChan, _ := vcf.GoReadToChan(vcfPath)
	for rec := range vcfChan {
		chrom := rec.Chr
		pos := rec.Pos
		chromMap, ok := divs[chrom]
		if !ok {
			continue
		}
		divList, ok := chromMap[pos]
		if !ok {
			continue
		}
		for _, d := range divList {
			anc := d.AncBase
			if requireSnp && len(rec.Ref) != 1 {
				continue
			}
			desiredAlleleIdx := -1
			matchedAlt := "."
			if anc == rec.Ref {
				desiredAlleleIdx = 0
			} else {
				for ai, alt := range rec.Alt {
					if strings.ToUpper(alt) == anc {
						desiredAlleleIdx = ai + 1
						matchedAlt = alt
						break
					}
				}
			}
			if desiredAlleleIdx == -1 {
				continue
			}
			ac, an := alleleCountsForVcf(rec, desiredAlleleIdx)
			if an == 0 {
				continue
			}
			af := float64(ac) / float64(an)
			if af > minFreq {
				vcfAltOut := matchedAlt
				if anc == rec.Ref {
					vcfAltOut = "."
				}
				start0 := d.Start0
				fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\n",
					d.Chrom, start0, rec.Id, rec.Ref, vcfAltOut, af, d.Name, d.DivBase, d.AncBase)
			}
		}
	}
}

func main() {
	var divFile string
	var vcfFile string
	var vcfDir string
	var outFile string
	var minFreq float64
	var requireSnp bool

	flag.StringVar(&divFile, "div", "", "Input div BED file (chrom start end name score strand div anc). Required.")
	flag.StringVar(&vcfFile, "vcf", "", "Single bgzipped VCF file to use (mutually exclusive with -vcfDir).")
	flag.StringVar(&vcfDir, "vcfDir", "", "Directory containing per-chrom VCFs (files like chr22.segSite.vcf.gz).")
	flag.StringVar(&outFile, "out", "", "Output file (default stdout).")
	flag.Float64Var(&minFreq, "minFreq", 0.01, "Minimum ancestral allele frequency threshold (default 0.01).")
	flag.BoolVar(&requireSnp, "snponly", true, "Only consider SNP records (default true).")
	flag.Parse()

	if divFile == "" {
		flag.Usage()
		log.Fatal("Error: -div is required.")
	}
	if vcfFile == "" && vcfDir == "" {
		flag.Usage()
		log.Fatal("Error: supply either -vcf or -vcfDir.")
	}
	if vcfFile != "" && vcfDir != "" {
		flag.Usage()
		log.Fatal("Error: supply only one of -vcf or -vcfDir.")
	}

	divs, err := readDivBed(divFile)
	if err != nil {
		log.Fatalf("failed reading div bed file: %v", err)
	}

	var outHandle *os.File
	if outFile == "" {
		outHandle = os.Stdout
	} else {
		outHandle, err = os.Create(outFile)
		if err != nil {
			log.Fatalf("could not create output file %s: %v", outFile, err)
		}
		defer outHandle.Close()
	}

	// header
	fmt.Fprintf(outHandle, "chr\tstart\tvcfId\tvcfRefBase\tvcfAltBase\tvcfFreqOfHaqerAnc\thaqerId\thaqerDivBase\thaqerAncBase\n")

	if vcfFile != "" {
		processUsingSingleVcf(divs, vcfFile, minFreq, requireSnp, outHandle)
	} else {
		processUsingVcfDir(divs, vcfDir, minFreq, requireSnp, outHandle)
	}
}
