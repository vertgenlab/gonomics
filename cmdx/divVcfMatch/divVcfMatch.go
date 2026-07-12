// Notes:
// CLI inputs below
// -div == files like /work/yl726/divsOutput.bed
// -vcfDir == /hpc/group/vertgenlab/raven/PrimateT2T_20way/selection_hs1/populationData/Filtered, files like chr22.segSite.vcf.gz
// -out == file like /work/yl726/divsVcfCheck.bed

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

// DivRecord represents one line from your divergence output file.
type DivRecord struct {
	Chrom string
	Start int // 0-based (as in your divergence file)
	End   int
	Div   string
	Anc   string
	Name  string
}

// map key type: chrom -> pos1 -> []DivRecord
type DivMap map[string]map[int][]DivRecord

func readDivFile(path string) (DivMap, error) {
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
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 6 {
			return nil, fmt.Errorf("line %d: expected >=6 columns (chrom start end div anc name), got %d: %q", lineNo, len(fields), line)
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
		div := fields[3]
		anc := strings.ToUpper(fields[4])
		name := fields[5]

		rec := DivRecord{
			Chrom: chrom,
			Start: start0,
			End:   end0,
			Div:   div,
			Anc:   anc,
			Name:  name,
		}
		pos1 := start0 + 1 // convert 0-based start to 1-based VCF POS

		if _, ok := ans[chrom]; !ok {
			ans[chrom] = make(map[int][]DivRecord)
		}
		ans[chrom][pos1] = append(ans[chrom][pos1], rec)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return ans, nil
}

// compute counts (AC, AN) of a desired allele index in a Vcf record.
// desiredAlleleIdx: 0 == REF, 1..n == ALT index+1
func alleleCountsForVcf(v vcf.Vcf, desiredAlleleIdx int) (ac int, an int) {
	ac = 0
	an = 0

	for _, s := range v.Samples {
		// prefer pre-parsed sample alleles if present
		if len(s.Alleles) > 0 {
			for _, a := range s.Alleles {
				if int(a) >= 0 { // -1 may indicate missing
					an++
					if int(a) == desiredAlleleIdx {
						ac++
					}
				}
			}
			continue
		}

		// fallback: try to parse GT from FormatData[0]
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
		// no genotype info for this sample; skip.
	}
	return ac, an
}

func processChrom(chrom string, chromMap map[int][]DivRecord, vcfDir string, minFreq float64, requireSnp bool, out *os.File) {
	vcfName := chrom + ".segSite.vcf.gz"
	vcfPath := filepath.Join(vcfDir, vcfName)

	if _, err := os.Stat(vcfPath); os.IsNotExist(err) {
		log.Printf("warning: VCF for %s not found at %s — skipping\n", chrom, vcfPath)
		return
	}

	vcfChan, _ := vcf.GoReadToChan(vcfPath)

	for rec := range vcfChan {
		pos := rec.Pos
		divRecs, ok := chromMap[pos]
		if !ok {
			continue
		}

		for _, d := range divRecs {
			anc := d.Anc
			if requireSnp && len(rec.Ref) != 1 {
				continue
			}

			desiredAlleleIdx := -1
			if anc == rec.Ref {
				desiredAlleleIdx = 0
			} else {
				for ai, alt := range rec.Alt {
					if anc == alt {
						desiredAlleleIdx = ai + 1
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
				fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%.6f\t%d\t%d\t%s\t%s\t%s\n",
					chrom, pos, d.Name, anc, af, ac, an, rec.Id, rec.Ref, strings.Join(rec.Alt, ","))
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

	flag.StringVar(&divFile, "div", "", "Divergence file (tab: chrom start end div anc name) (required)")
	flag.StringVar(&vcfFile, "vcf", "", "Single bgzipped VCF file to use (mutually exclusive with -vcfDir)")
	flag.StringVar(&vcfDir, "vcfDir", "", "Directory containing per-chrom VCFs (files like chr22.segSite.vcf.gz)")
	flag.StringVar(&outFile, "out", "", "Output file (default stdout)")
	flag.Float64Var(&minFreq, "minFreq", 0.01, "Minimum ancestral allele frequency threshold (default 0.01 for 1%)")
	flag.BoolVar(&requireSnp, "snponly", true, "Only consider SNP records (default true)")
	flag.Parse()

	if divFile == "" {
		flag.Usage()
		log.Fatal("Error: -div is required")
	}
	if vcfFile == "" && vcfDir == "" {
		flag.Usage()
		log.Fatal("Error: supply either -vcf (single file) or -vcfDir (directory of per-chrom VCFs)")
	}
	if vcfFile != "" && vcfDir != "" {
		flag.Usage()
		log.Fatal("Error: supply only one of -vcf or -vcfDir")
	}

	divMap, err := readDivFile(divFile)
	if err != nil {
		log.Fatalf("failed reading divergence file: %v", err)
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
	fmt.Fprintf(outHandle, "chrom\tpos\tname\tanc\tAF\tAC\tAN\tvcfId\tREF\tALT\n")

	// If user supplied a single VCF file, process only the chromosomes that map to that file.
	if vcfFile != "" {
		vcfPath := vcfFile
		if _, err := os.Stat(vcfPath); os.IsNotExist(err) {
			log.Fatalf("vcf file not found: %s", vcfPath)
		}
		vcfChan, _ := vcf.GoReadToChan(vcfPath)
		// build a set of positions we care about from all chroms (but the VCF will have a CHROM field)
		// simpler: as we stream, check if rec.Chr is in divMap
		for rec := range vcfChan {
			chrom := rec.Chr
			if chromMap, ok := divMap[chrom]; ok {
				// reuse the chrom processing logic by creating a tiny temp map for this rec.Pos
				// but simpler: process the single rec inline:
				divRecs, ok2 := chromMap[rec.Pos]
				if !ok2 {
					continue
				}
				for _, d := range divRecs {
					anc := d.Anc
					if requireSnp && len(rec.Ref) != 1 {
						continue
					}
					desiredAlleleIdx := -1
					if anc == rec.Ref {
						desiredAlleleIdx = 0
					} else {
						for ai, alt := range rec.Alt {
							if anc == alt {
								desiredAlleleIdx = ai + 1
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
						fmt.Fprintf(outHandle, "%s\t%d\t%s\t%s\t%.6f\t%d\t%d\t%s\t%s\t%s\n",
							chrom, rec.Pos, d.Name, anc, af, ac, an, rec.Id, rec.Ref, strings.Join(rec.Alt, ","))
					}
				}
			}
		}
	} else {
		// vcfDir mode: iterate chromosomes in divMap
		for chrom, chromMap := range divMap {
			processChrom(chrom, chromMap, vcfDir, minFreq, requireSnp, outHandle)
		}
	}
}
