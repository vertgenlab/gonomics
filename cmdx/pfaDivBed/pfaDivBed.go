// Notes:
// Take a file of bed regions, e.g. chr1 100 200
// For each bed region, read the appropriate pFa file, e.g. chr1.pFa, extract subsequences e.g. positions 100-200, then scan 2 species in the []pFa, e.g. firstQuery, secondQuery, for divergent positions.
// For each divergent position, report chrom, start, end, in reference pFa coordinates, as well as most likely bases in firstQuery and secondQuery.
// Later, I will compare each bed region's divergent positions with population data vcf.
// did not add safety check code
// for my purposes, the pfa is multi-pfa, all from 1 chr, but many species
// for my purposes, firstQueryName == HUMANanc (Div), secondQueryName == hcaT2T (Anc)
// hard-coded filename .pFa.gz
// DivBed Name will be HAQER name/ID
// CLI inputs below
// bedFile == /hpc/group/vertgenlab/raven/PrimateT2T_20way/aqers/T2T/niceName/human.bed, Consensus HAQER bed, 1596 lines, column 4 is HAQER name/ID
// pfaDir == /hpc/group/vertgenlab/raven/PrimateT2T_20way/work/yl726/PrimateT2T_20way/haqerTmp.humanT2T, files like chr22.pFa.gz

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
	"path/filepath"
	"strings"
)

// Example of integrating into main: add flags for bed and pfa directory and call processRegionsFromPerChromPfa
type Settings struct {
	BedFile                 string
	PfaDir                  string
	OutFile                 string
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

// processRegionsFromPerChromPfa processes a BED file, for each region opens a chromosome-specific pFa file (chrom + ".pFa.gz"),
// extracts the window [ChromStart:ChromEnd) from the pFasta alignment, selects reference/fistQuery/secondQuery sequences,
// and accumulates divergence calls using pfaDivBed.
// - pfaDir is the dir or prefix where per-chrom pFa files live (can be "" if bed regions point to full path).
// - firstQueryName/secondQueryName are the names of the two species/records inside the pFa file. If empty, fall back to records[1] and records[2].
// - baseDotThreshold is the threshold for calling a substitution.
func processRegionsFromPerChromPfa(bedFile string, pfaDir string, firstQueryName string, secondQueryName string, baseDotThreshold float64) (allDivs []DivBed) {
	regions := bed.Read(bedFile)

	// for each bed region
	for _, reg := range regions {

		// Step 1: get pFa file
		// build per-chrom pFa path; allow reg.Chrom to already contain path if pfaDir empty or user passed full path in bed.
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

		// Step 3: from 3 pFa, slice out region [ChromStart:ChromEnd)
		if reg.ChromStart < 0 || reg.ChromEnd > len(reference) {
			log.Fatalf("Error: region %s:%d-%d out of range for %s (len=%d)\n", reg.Chrom, reg.ChromStart, reg.ChromEnd, pfaName, len(reference))
		}

		refSub := reference[reg.ChromStart:reg.ChromEnd]
		firstSub := firstQuery[reg.ChromStart:reg.ChromEnd]
		secondSub := secondQuery[reg.ChromStart:reg.ChromEnd]

		// Step 4: for sliced out 3 pFa [ChromStart:ChromEnd), get divs
		divs := pfaDivBed(refSub, firstSub, secondSub, baseDotThreshold, reg.Chrom, reg.ChromStart, reg.Name)

		// Step 5: append divs for this bed region to allDivs
		allDivs = append(allDivs, divs...)
	}

	return allDivs
}

// Step 4 helper function
// pfaDivBed scans aligned slices for substitutions (dot-product based) and returns DivBed entries.
// reference, firstQuery, secondQuery should be alignment slices corresponding to the same alignment window.
// baseDotToSubstThreshold is the threshold above which a position is called a substitution.
// chrom is the chromosome name for BED output.
// refStartOffset is the reference coordinate offset (e.g., region.ChromStart) to convert local indices to absolute reference coordinates.
func pfaDivBed(referenceSub []pDna.Float32Base, firstQuerySub []pDna.Float32Base, secondQuerySub []pDna.Float32Base, baseDotToSubstThreshold float64, chrom string, refStartOffset int, name string) (divAnswers []DivBed) {
	if len(referenceSub) != len(firstQuerySub) || len(referenceSub) != len(secondQuerySub) {
		log.Fatalf("Error: sequences passed to pfaDivBed must be equal length.")
	}

	var start, end int
	var div, anc string
	var baseDot float64
	lastRefStart := 0
	lastAlnStart := 0

	// copied code from pfaFindFast "8. numSubst: is this a substitution?"
	for alnIdx := 0; alnIdx < len(firstQuerySub); alnIdx++ {

		// skip positions where either sequence has a gap
		if !pDna.IsGap(firstQuerySub[alnIdx]) && !pDna.IsGap(secondQuerySub[alnIdx]) { // do not calculate over gaps
			baseDot = pDna.DotSubstProb(firstQuerySub[alnIdx], secondQuerySub[alnIdx]) // for non-gap position, calculate substitution probability from dot product method
			if baseDot >= baseDotToSubstThreshold {                                    // if the substitution probability >= threshold
				// then this position is a substitution

				// map alignment index to reference position (absolute)
				// use PAlnPosToRefPosCounterSeq and pass refStart = refStartOffset + lastRefStart
				start = pFasta.PAlnPosToRefPosCounterSeq(referenceSub, alnIdx, lastRefStart, lastAlnStart) + refStartOffset
				// update counters so mapping stays efficient
				lastAlnStart = alnIdx
				lastRefStart = start - refStartOffset

				end = start + 1 // BED end is exclusive; will add +0 here; when printing add +1 if you want 1-based end

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

// Step 6: write to output file
// writeDivBedsToFile writes Divergent positions in BED-like tabular format:
// chrom  start  end  divBase  ancBase	name
// NOTE: this writes end as exclusive. If you want a 1-based half-open BED, you can +1 end here.
func writeDivBedsToFile(fileName string, records []DivBed) {
	var err error

	//create file
	file := fileio.EasyCreate(fileName)

	//loop through structs and write each struct to a line
	for _, r := range records {
		line := fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s", r.Chrom, r.StartPos, r.EndPos, r.Div, r.Anc, r.Name)
		fileio.WriteToFileHandle(file, line)
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"pfaDivBed \n" +
			"Usage:\n" +
			" pfaFindFast input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	pfaDir := flag.String("pfaDir", ".", "Directory where per-chrom pFa files live (each named CHROM.pfa).")
	firstQ := flag.String("firstQueryName", "", "Name of first query in pFa.")
	secondQ := flag.String("secondQueryName", "", "Name of second query in pFa.")
	baseDot := flag.Float64("baseDotToSubstThreshold", 0.8, "Threshold for dot-product substitution call.")
	flag.Parse()

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	s := Settings{
		BedFile:                 flag.Arg(0),
		PfaDir:                  *pfaDir,
		OutFile:                 flag.Arg(1),
		FirstQueryName:          *firstQ,
		SecondQueryName:         *secondQ,
		BaseDotToSubstThreshold: *baseDot,
	}

	divs := processRegionsFromPerChromPfa(s.BedFile, s.PfaDir, s.FirstQueryName, s.SecondQueryName, s.BaseDotToSubstThreshold)
	writeDivBedsToFile(s.OutFile, divs)
}
