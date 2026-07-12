// Notes:
// CLI inputs below
// bedFile == /hpc/group/vertgenlab/raven/PrimateT2T_20way/aqers/T2T/niceName/human.bed, Consensus HAQER bed, 1596 lines, column 4 is HAQER name/ID
// pfaDir == /hpc/group/vertgenlab/raven/PrimateT2T_20way/work/yl726/PrimateT2T_20way/haqerTmp.humanT2T, files like chr22.pFa.gz
// outDir == /work/yl726/pfaDivBedOutputs

// More Notes:
// Take a file of bed regions, e.g. chr1 100 200
// For each bed region, read the appropriate pFa file, e.g. chr1.pFa, extract subsequences e.g. positions 100-200, then scan 2 species in the []pFa, e.g. firstQuery, secondQuery, for divergent positions.
// For each divergent position, report chrom, start, end, in reference pFa coordinates, as well as most likely bases in firstQuery and secondQuery.
// Later, I will compare each bed region's divergent positions with population data vcf.
// did not add safety check code
// for my purposes, the pfa is multi-pfa, all from 1 chr, but many species
// for my purposes, firstQueryName == HUMANanc (Div), secondQueryName == hcaT2T (Anc)
// hard-coded filename .pFa.gz
// DivBed Name will be HAQER name/ID

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"path/filepath"
	"strings"
)

// Example of integrating into main: add flags for bed and pfa directory and call processRegionsFromPerChromPfa
type Settings struct {
	BedFile                 string
	PfaDir                  string
	OutDir                  string
	FirstQueryName          string
	SecondQueryName         string
	BaseDotToSubstThreshold float64
}

// DivBed: output struct for divergent positions
type DivBed struct {
	Chrom    string
	StartPos int
	EndPos   int
	Div      string // most likely base in firstQuery (divergent)
	Anc      string // most likely base in secondQuery (ancestral)
	Name     string // record name
}

// convert a reference position (0-based) to an alignment position (0-based) in a pFasta sequence.
// works like fasta.RefPosToAlnPosCounter / RefPosToAlnPos but for []pDna.Float32Base
func pRefPosToAlnPosCounterSeq(record []pDna.Float32Base, refPos int, refStart int, alnStart int, allowBed bool) int {
	if refStart > refPos {
		log.Fatalf("refStart > RefPos")
	}
	if alnStart >= len(record) {
		log.Fatalf("Ran out of alignment.")
	}

	// behaviour mirrors fasta.refPosToAlnPosCounterExposed (allowBed controls end-exclusive safety)
	if !allowBed {
		for t := alnStart; refStart < refPos; alnStart++ {
			t++
			if t == len(record) {
				log.Fatalf("Ran out of chromosome.")
			} else if !pDna.IsGap(record[t]) {
				refStart++
			}
		}
	} else {
		// allowBed: permit stepping one past the end for bed end coordinates
		for t := alnStart; refStart < refPos; alnStart++ {
			t++
			if t > len(record) {
				log.Fatalf("Ran out of chromosome.")
			} else if t == len(record) {
				// push alnStart one forward for the bed-end convention then break
				alnStart++
				break
			} else if !pDna.IsGap(record[t]) {
				refStart++
			}
		}
	}
	return alnStart
}

func pRefPosToAlnPosSeq(record []pDna.Float32Base, refPos int) int {
	return pRefPosToAlnPosCounterSeq(record, refPos, 0, 0, false)
}

func pRefPosToAlnPosBedSeq(record []pDna.Float32Base, refPos int) int {
	return pRefPosToAlnPosCounterSeq(record, refPos, 0, 0, true)
}

// processRegionsFromPerChromPfa processes a BED file, for each region opens a chromosome-specific pFa file (chrom + ".pFa.gz"),
// extracts the window [ChromStart:ChromEnd) from the pFasta alignment, selects reference/fistQuery/secondQuery sequences,
// and writes a per-region bed file named <Name>.bed into outDir.
func processRegionsFromPerChromPfa(bedFile string, pfaDir string, outDir string, firstQueryName string, secondQueryName string, baseDotThreshold float64) {
	regions := bed.Read(bedFile)

	// for each bed region
	for _, reg := range regions {

		// Step 1: get pFa file
		pfaName := reg.Chrom + ".pFa.gz"
		if pfaDir != "" {
			pfaName = filepath.Join(pfaDir, pfaName)
		}

		records := pFasta.Read(pfaName)

		// Step 2: from pFa file, get firstQuery, secondQuery, reference pFa
		if len(records) < 2 {
			log.Fatalf("Error: pFasta %s must contain at least 2 records (reference + queries).\n", pfaName)
		}
		recordsMap := pFasta.ToMap(records)

		var reference, firstQuery, secondQuery []pDna.Float32Base
		reference = records[0].Seq // convention: first record is reference

		// load firstQuery
		if strings.TrimSpace(firstQueryName) == "" {
			// fallback: use record 1 (if exists)
			firstQuery = records[1].Seq
		} else {
			var ok bool
			firstQuery, ok = recordsMap[firstQueryName]
			if !ok {
				log.Fatalf("Error: firstQueryName '%s' not found in %s\n", firstQueryName, pfaName)
			}
		}

		// load secondQuery
		if strings.TrimSpace(secondQueryName) == "" {
			// fallback: if available use record 2 else error
			if len(records) < 3 {
				log.Fatalf("Error: secondQueryName not specified and pFasta %s contains fewer than 3 records.\n", pfaName)
			}
			secondQuery = records[2].Seq
		} else {
			var ok bool
			secondQuery, ok = recordsMap[secondQueryName]
			if !ok {
				log.Fatalf("Error: secondQueryName '%s' not found in %s\n", secondQueryName, pfaName)
			}
		}

		// ensure alignment lengths equal
		if !(len(reference) == len(firstQuery) && len(reference) == len(secondQuery)) {
			log.Fatalf("Error: alignment lengths differ in %s\n", pfaName)
		}

		// Step 3: convert reference coords to alignment coords, then slice by aln indices
		// Ensure supplied ref coordinates are in [0, referenceRefLength]
		refLength := 0
		// compute reference length (number of non-gap positions in the reference pFasta)
		for i := 0; i < len(reference); i++ {
			if !pDna.IsGap(reference[i]) {
				refLength++
			}
		}
		if reg.ChromStart < 0 || reg.ChromEnd > refLength {
			log.Fatalf("Error: region %s:%d-%d out of reference range [0,%d) for %s\n", reg.Chrom, reg.ChromStart, reg.ChromEnd, refLength, pfaName)
		}

		// map ref positions to alignment indices
		alnStart := pRefPosToAlnPosSeq(reference, reg.ChromStart)
		alnEnd := pRefPosToAlnPosBedSeq(reference, reg.ChromEnd) // allow bed-end behaviour

		if alnStart < 0 || alnEnd > len(reference) || alnStart >= alnEnd {
			log.Fatalf("Error mapping ref->aln for region %s:%d-%d -> aln %d-%d in %s\n", reg.Chrom, reg.ChromStart, reg.ChromEnd, alnStart, alnEnd, pfaName)
		}

		refSub := reference[alnStart:alnEnd]
		firstSub := firstQuery[alnStart:alnEnd]
		secondSub := secondQuery[alnStart:alnEnd]

		// Step 4: for sliced out 3 pFa [ChromStart:ChromEnd), get divs
		divs := pfaDivBed(refSub, firstSub, secondSub, baseDotThreshold, reg.Chrom, reg.ChromStart, reg.Name)

		// Step 5: write per-region bed file
		// sanitize name for filename
		safeName := sanitizeName(reg.Name)
		outPath := filepath.Join(outDir, safeName+".bed")
		writeDivBedsToFile(outPath, divs)
	}
}

// pfaDivBed scans aligned slices for substitutions and returns DivBed entries.
func pfaDivBed(referenceSub []pDna.Float32Base, firstQuerySub []pDna.Float32Base, secondQuerySub []pDna.Float32Base, baseDotToSubstThreshold float64, chrom string, refStartOffset int, name string) (divAnswers []DivBed) {
	if len(referenceSub) != len(firstQuerySub) || len(referenceSub) != len(secondQuerySub) {
		log.Fatalf("Error: sequences passed to pfaDivBed must be equal length.")
	}

	var start, end int
	var div, anc string
	var baseDot float64
	lastRefStart := 0
	lastAlnStart := 0

	for alnIdx := 0; alnIdx < len(firstQuerySub); alnIdx++ {
		// skip positions where either sequence has a gap
		if !pDna.IsGap(firstQuerySub[alnIdx]) && !pDna.IsGap(secondQuerySub[alnIdx]) {
			baseDot = pDna.DotSubstProb(firstQuerySub[alnIdx], secondQuerySub[alnIdx])
			if baseDot >= baseDotToSubstThreshold {
				// map alignment index to reference position (slice-relative), then add offset
				startLocal := pFasta.PAlnPosToRefPosCounterSeq(referenceSub, alnIdx, lastRefStart, lastAlnStart)
				start = startLocal + refStartOffset
				// update counters (keep lastRefStart relative to slice)
				lastAlnStart = alnIdx
				lastRefStart = startLocal

				end = start + 1

				div = dna.BaseToString(pDna.MostLikelyBase(firstQuerySub[alnIdx]))
				anc = dna.BaseToString(pDna.MostLikelyBase(secondQuerySub[alnIdx]))

				divAnswers = append(divAnswers, DivBed{
					Chrom:    chrom,
					StartPos: start,
					EndPos:   end,
					Div:      div,
					Anc:      anc,
					Name:     name,
				})
			}
		}
	}
	return divAnswers
}

// writeDivBedsToFile writes Divergent positions in BED-like tabular format:
// chrom  start  end  divBase  ancBase  name
func writeDivBedsToFile(fileName string, records []DivBed) {
	file := fileio.EasyCreate(fileName)
	defer func() {
		err := file.Close()
		exception.PanicOnErr(err)
	}()

	for _, r := range records {
		line := fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\n", r.Chrom, r.StartPos, r.EndPos, r.Div, r.Anc, r.Name)
		fileio.WriteToFileHandle(file, line)
	}
	// If records is empty, we still created the file (empty). Change behavior here if you want to skip empty files.
}

// sanitizeName replaces spaces and slashes with safe characters for filenames.
func sanitizeName(name string) string {
	s := strings.TrimSpace(name)
	s = strings.ReplaceAll(s, " ", "_")
	s = strings.ReplaceAll(s, "/", "_")
	s = strings.ReplaceAll(s, "\\", "_")
	return s
}

func usage() {
	fmt.Print(
		"pfaDivBed \n" +
			"Usage:\n" +
			" pfaDivBed input.bed outDir\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	pfaDir := flag.String("pfaDir", ".", "Directory where per-chrom pFa files live (each named CHR.pFa.gz).")
	outDir := flag.String("outDir", ".", "Directory to write per-region bed files (default current directory).")
	firstQ := flag.String("firstQueryName", "", "Name of first query in pFa.")
	secondQ := flag.String("secondQueryName", "", "Name of second query in pFa.")
	baseDot := flag.Float64("baseDotToSubstThreshold", 0.8, "Threshold for dot-product substitution call.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: expecting bed input file as positional argument.")
	}
	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: expecting at least 1 positional argument.")
	}

	s := Settings{
		BedFile:                 flag.Arg(0),
		PfaDir:                  *pfaDir,
		OutDir:                  *outDir,
		FirstQueryName:          *firstQ,
		SecondQueryName:         *secondQ,
		BaseDotToSubstThreshold: *baseDot,
	}

	// Create outDir
	var err error
	if mkerr := os.MkdirAll(s.OutDir, 0755); mkerr != nil {
		log.Fatalf("could not create outDir %s: %v (fileio.CreateDirsIfNeeded err: %v)", s.OutDir, mkerr, err)
	}

	processRegionsFromPerChromPfa(s.BedFile, s.PfaDir, s.OutDir, s.FirstQueryName, s.SecondQueryName, s.BaseDotToSubstThreshold)
}
