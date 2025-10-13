package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"path/filepath"
	"strings"
)

type Settings struct {
	FaDir   string //input fasta directory with all chromosomes
	Bed     string //input bed file
	BaseOne string //base one in pair
	BaseTwo string //base two in pair
	Outfile string //output file name
	Compare bool   //-compare mode for 2 aligned sequences
}

type PairOfBaseCounts struct {
	Chrom        string
	Name         string //region name from bed file
	Gain         int
	Loss         int
	Cons         int
	WeakToStrong int //subs
	StrongToWeak int //subs
}

// -compare mode: check if constant pair
func isCons(firstSeq1, firstSeq2, secondSeq1, secondSeq2, b1, b2 dna.Base) bool {
	if (firstSeq1 == b1 && firstSeq2 == b2) && (secondSeq1 == b1 && secondSeq2 == b2) {
		return true
	} else {
		return false
	}
}

// -compare mode: check if gained pair in sequence 1
func isGain(firstSeq1, firstSeq2, secondSeq1, secondSeq2, b1, b2 dna.Base) bool {
	if (firstSeq1 == b1 && firstSeq2 == b2) && !(secondSeq1 == b1 && secondSeq2 == b2) {
		return true
	} else {
		return false
	}
}

// -compare mode: check if lost pair in sequence 1
func isLoss(firstSeq1, firstSeq2, secondSeq1, secondSeq2, b1, b2 dna.Base) bool {
	if !(firstSeq1 == b1 && firstSeq2 == b2) && (secondSeq1 == b1 && secondSeq2 == b2) {
		return true
	} else {
		return false
	}
}

// -compare mode: track strong -> weak substitutions (C||G -> A||T)
func isStrongToWeak(firstSeq1, secondSeq1 dna.Base) bool {
	if (firstSeq1 == dna.A || firstSeq1 == dna.T) && (secondSeq1 == dna.C || secondSeq1 == dna.G) {
		return true
	} else {
		return false
	}
}

// -compare mode: track weak -> strong substitutions (A||T -> C||G)
func isWeakToStrong(firstSeq1, secondSeq1 dna.Base) bool {
	if (firstSeq1 == dna.C || firstSeq1 == dna.G) && (secondSeq1 == dna.A || secondSeq1 == dna.T) {
		return true
	} else {
		return false
	}
}

// as you're scanning the sequence, nextBase() recognizes when base 2 in the pair is a gap, traverses the gaps until a base is found, and uses that as "base 2" instead
func nextBase(region []dna.Base, currPos int) (nextBase dna.Base) {
	for i := currPos; i < len(region); i++ {
		if dna.DefineBase(region[i]) {
			return region[i]
		}
	}
	return dna.Gap
}

// find correct chromosome fasta for the BED regions being analyzed
func findChromInFile(fileDir string, chrom string) (string, error) {
	//09-29-25: force naming convention: chr1.fa, chr2.fa, etc. so chr name can be extracted easily
	files, err := os.ReadDir(fileDir)
	if err != nil {
		return "", err
	}
	expectedName := chrom + ".fa"
	for _, f := range files {
		fileName := f.Name()
		if fileName == expectedName || fileName == expectedName+".gz" {
			filePath := filepath.Join(fileDir, fileName)
			return filePath, nil
		}
	}
	return "", fmt.Errorf("%s fasta not found in %s", chrom, fileDir)
}

// RefPosToAlnPosBed loops existing functions for converting reference positions to alignment positions
// (RefPosToAlnPos and RefPostoAlnPosCount), to convert the coordinates of an entire bed file.
// It takes a bed file and a multiFa as inputs and outputs the bed regions with alignment coordinates.

func RefPosToAlnPosBed(input []bed.Bed, fastaFile string) (output []bed.Bed) {
	var regName string

	seqs := fasta.Read(fastaFile)
	refSeq := seqs[0]
	firstRegion := input[0]
	chrom := firstRegion.Chrom

	if firstRegion.Name != "" {
		regName = firstRegion.Name
	} else {
		log.Fatalf("Error: each BED region must have a name in column 4")
	}

	//convert first region's start coordinate from ref to alignment pos
	startAlnPos := fasta.RefPosToAlnPos(refSeq, firstRegion.ChromStart)

	//convert first region's end coordinate from ref to alignment pos
	endAlnPos := fasta.RefPosToAlnPos(refSeq, firstRegion.ChromEnd)

	//add this region (chr, alignment start, alignment end, name) to slice of bed structs
	output = append(output, bed.Bed{Chrom: chrom, ChromStart: startAlnPos, ChromEnd: endAlnPos, Name: regName, FieldsInitialized: 4})

	//use the end coordinates of firstRegion (reference and alignment) for next function
	lastRefPos := firstRegion.ChromEnd
	lastAlnPos := endAlnPos

	//starting from second region, loop through regions in input bed
	for i := 1; i < len(input); i++ {
		region := input[i]

		if region.Name != "" {
			regName = region.Name
		} else {
			log.Fatalf("Error: each BED region must have a name in column 4")
		}

		//use the last reference/alignment position pair to convert the current region's ref start coordinate to an alignment start coordinate
		startAlnPos = fasta.RefPosToAlnPosCounter(refSeq, region.ChromStart, lastRefPos, lastAlnPos)

		//use the same strategy for the current region's end coordinate
		endAlnPos = fasta.RefPosToAlnPosCounter(refSeq, region.ChromEnd, lastRefPos, lastAlnPos)

		//add this information to output
		output = append(output, bed.Bed{Chrom: chrom, ChromStart: startAlnPos, ChromEnd: endAlnPos, Name: regName, FieldsInitialized: 4})

		//reset the ref/alignment position pair to that of the current region, so that it can be used for the next region
		lastRefPos = region.ChromEnd
		lastAlnPos = endAlnPos
	}
	return output
}

// count pairs of user-provided bases (such as CG or GC) in a SINGLE sequence fasta (default mode)
// inputs: fastaFile for the correct chromosome returned by findChromInFile, original reference bed file, base 1 and base 2 to count as a pair
// output: counts of that pair of bases for each genomic region
func countPairOfBase(fastaFile string, refBed []bed.Bed, b1, b2 dna.Base) (regionPairs []bed.Bed, err error) {
	chr := fasta.Read(fastaFile)

	//check that there are not more than 1 fasta sequence
	if len(chr) != 1 {
		log.Fatalf("Error: expecting exactly one record in fasta file, but got %d. If you want to compare 2 sequences, use --compare mode.\n", len(chr))
	}

	//grab seq
	seq := chr[0].Seq

	//if sequence is empty
	if len(seq) == 0 {
		log.Fatalf("Error: fasta sequence is empty.\n")
	}

	//go through BED regions
	for _, region := range refBed {
		var f1, f2 dna.Base
		//get the following information
		chrom := region.Chrom
		start := region.ChromStart
		end := region.ChromEnd
		name := region.Name

		//Find the current region in the fasta sequence:
		regSeq := seq[start:end]

		//initialize number of paired sites to 0
		pairOfBaseCount := 0

		//for length of entire sequence defined by reference start/end coordinates
		for i := 0; i < len(regSeq)-1; i++ {
			f1, f2 = regSeq[i], regSeq[i+1]
			if f1 == b1 && f2 == b2 {
				pairOfBaseCount++
			}
		}
		//output
		regionPairs = append(regionPairs, bed.Bed{Chrom: chrom, ChromStart: start, ChromEnd: end, Name: name, Score: pairOfBaseCount, FieldsInitialized: 5})
	}
	return regionPairs, nil
}

// -compare mode counting gains, losses, and constant pairs of bases provided by user
// inputs: multiFasta for the correct chromosome returned by findChromInFile, alignment-converted BED returned by RefPosToAlnPosBed, base 1 and base 2 to count as a pair
func comparePairOfBaseCount(fastaFile string, alnBed []bed.Bed, b1, b2 dna.Base) (regionPairs []PairOfBaseCounts, err error) {
	//get fasta entries
	seqs := fasta.Read(fastaFile)

	//if not 2 entries
	if len(seqs) != 2 {
		log.Fatalf("Error: expecting exactly two records in fasta file, but got %d. --compare mode compares exactly 2 aligned sequences.\n", len(seqs))
	}

	//define sequence 1 and sequence 2 in the multiFa
	firstSeq := seqs[0].Seq
	secondSeq := seqs[1].Seq

	//for each region in the alignment-converted bed file
	for _, region := range alnBed {
		//get the following information
		chrom := region.Chrom
		start := region.ChromStart
		end := region.ChromEnd
		name := region.Name

		//Find the current region in the fasta sequences:
		//(Here, we assume that both genome 1 and genome 2 are unmodified
		//from the original multiple alignment, so the same alignment coordinates apply to both)
		firstRegseq := firstSeq[start:end]
		secondRegseq := secondSeq[start:end]

		//initialize number of each category of paired sites to 0
		gainCount := 0
		lossCount := 0
		consCount := 0

		//initialize each category of substitutions to 0
		weakToStrongCount := 0
		strongToWeakCount := 0

		//for length of entire sequence defined by alignment start/end coordinates
		for i := 0; i < len(firstRegseq)-1; i++ {
			//f1, f2 and s1, s2 are calls to nextBase()
			f1, s1 := firstRegseq[i], secondRegseq[i]
			f2 := nextBase(firstRegseq, i+1)
			s2 := nextBase(secondRegseq, i+1)

			//check is constant, gained, or lost pair of bases
			switch {
			case isCons(f1, f2, s1, s2, b1, b2):
				consCount++
			case isGain(f1, f2, s1, s2, b1, b2):
				gainCount++
			case isLoss(f1, f2, s1, s2, b1, b2):
				lossCount++
			}

			//check for weak->strong or strong->weak substitutions
			if isWeakToStrong(f1, s1) {
				weakToStrongCount++
			}
			if isStrongToWeak(f1, s1) {
				strongToWeakCount++
			}

			//evaluate last base in the sequence for a weak->strong or strong->weak substitution
			if i == len(firstRegseq)-2 {
				if isWeakToStrong(f2, s2) {
					weakToStrongCount++
				}

				if isStrongToWeak(f2, s2) {
					strongToWeakCount++
				}
			}
		}
		//output
		regionPairs = append(regionPairs, PairOfBaseCounts{Chrom: chrom, Name: name, Gain: gainCount, Loss: lossCount, Cons: consCount, WeakToStrong: weakToStrongCount, StrongToWeak: strongToWeakCount})
	}
	return regionPairs, nil
}

// perform above^ functions to count pairs of bases in each region per chromosome
// this allows us to split up the BED file by chromosome internally, and only open up the corresponding chromosome Fasta file once
func countByChrom(s Settings) {
	var counts []bed.Bed
	var countsComp []PairOfBaseCounts
	var regCounts bed.Bed
	var regCountsComp PairOfBaseCounts

	//convert user-provided bases to dna.Base type
	baseOne := dna.StringToBase(strings.TrimSpace(s.BaseOne))
	baseTwo := dna.StringToBase(strings.TrimSpace(s.BaseTwo))

	bedBychrom := make(map[string][]bed.Bed)
	bedFile := s.Bed
	//preserve the original order of bed regions in the file
	bedRegions := bed.Read(bedFile)

	//make map of BED regions by chromosome, to process one chromosome at a time
	for _, region := range bedRegions {
		bedBychrom[region.Chrom] = append(bedBychrom[region.Chrom], region)
	}

	file := fileio.EasyCreate(s.Outfile)

	//write file header
	if !s.Compare {
		fileio.WriteToFileHandle(file, "Chrom\tStart\tEnd\tName\tPairOfBaseCount")
	} else {
		fileio.WriteToFileHandle(file, "Chrom\tStart\tEnd\tName\tGain\tLoss\tCons\tWeakToStrongSubs\tStrongToWeakSubs")
	}

	//fast lookup table, PER CHROMOSOME, to match paired base counts to the region's name in comp mode
	chromCountsMap := make(map[string]map[string]bed.Bed)
	//fast lookup table, PER CHROMOSOME, to match paired base counts to the region's name in comp mode
	chromCountsCompMap := make(map[string]map[string]PairOfBaseCounts)

	//for each chromosome
	for chrom, regions := range bedBychrom {
		//open the corresponding fasta
		fastaFile, err := findChromInFile(s.FaDir, chrom)
		exception.PanicOnErr(err)

		//if default mode
		if !s.Compare {
			//get counts per region
			counts, err = countPairOfBase(fastaFile, regions, baseOne, baseTwo)
			exception.PanicOnErr(err)

			//make map for these counts keyed by name for easy lookup
			cpgMap := make(map[string]bed.Bed)

			for _, r := range counts {
				cpgMap[r.Name] = r
			}
			//put those in a larger map by chrom
			chromCountsMap[chrom] = cpgMap
		} else {
			//if in -compare mode
			//convert BED to alignment coordinates
			alnBed := RefPosToAlnPosBed(regions, fastaFile)

			//get counts
			countsComp, err = comparePairOfBaseCount(fastaFile, alnBed, baseOne, baseTwo)
			exception.PanicOnErr(err)

			//make maps as before
			cpgMap := make(map[string]PairOfBaseCounts)

			for _, c := range countsComp {
				cpgMap[c.Name] = c
			}
			chromCountsCompMap[chrom] = cpgMap
		}
	}

	// output regions in original BED order by looping through the regions from the original BED file
	for _, region := range bedRegions {

		//if default mode
		if !s.Compare {
			//access all of the region counts for the region's chromosome
			countsMap, ok := chromCountsMap[region.Chrom]
			if !ok {
				log.Fatalf("No counts found for chromosome %s", region.Chrom)
			}

			//then pull out the counts for JUST this region
			regCounts, ok = countsMap[region.Name]
			if !ok {
				log.Fatalf("Counts not found for region %s", region.Name)
			}
		} else {
			//if in -compare mode
			//access all of the region counts for the region's chromosome
			countsCompMap, ok := chromCountsCompMap[region.Chrom]
			if !ok {
				log.Fatalf("No counts found for chromosome %s", region.Chrom)
			}

			//then pull out the counts for JUST this region
			//need to be separate from default mode because they are different structs
			regCountsComp, ok = countsCompMap[region.Name]
			if !ok {
				log.Fatalf("Counts not found for region %s", region.Name)
			}
		}

		//if default mode
		if !s.Compare {
			//print info from all data structs related to that region:
			//region.Chrom, region.ChromStart, etc. all come from the BED entry in the original BED file
			//regCounts.Score are the paired base coutns that were just pulled earlier^ based on the region's name (so we know they're correct)
			line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d",
				region.Chrom, region.ChromStart, region.ChromEnd, region.Name, regCounts.Score)
			fileio.WriteToFileHandle(file, line)
		} else {
			//if -compare mode
			//print info from all data structs related to that region:
			//region.Chrom, region.ChromStart, etc. all come from the BED entry in the original BED file
			//regCountsComp.Gain, Loss, etc. are the paired base counts/substitution counts that were just pulled earlier^ based on the region's name (so we know they're correct)
			line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d",
				region.Chrom, region.ChromStart, region.ChromEnd, region.Name,
				regCountsComp.Gain, regCountsComp.Loss, regCountsComp.Cons, regCountsComp.WeakToStrong, regCountsComp.StrongToWeak)
			fileio.WriteToFileHandle(file, line)
		}
	}

	//close file
	err := file.Close()
	exception.PanicOnErr(err)

	fmt.Println("Pair counts found and written to", s.Outfile)
}

func usage() {
	fmt.Print("countPairOfBases [options] fastaDir bedfileName baseOne baseTwo outfileName\n" +
		"\tCounts pairs of bases (ex: CG) within BED regions.\n" +
		"\tARGS: \n" +
		"\tfastaDir - path to directory containing multiFas for all chromosomes to analyze. Files must be named: {chromosome name}.fa (files can be gzipped).\n" +
		"\tbedfileName - name of BED file containing all genomic regions in which to count CpG changes. All regions are required to have a name in column 4.\n" +
		"\t**Note: please ensure that all chromosomes in the BED file have corresponding fastas in the provided fasta directory.**\n" +
		"\tbaseOne - first base of pair to find (ex: 'C')\n" +
		"\tbaseTwo - second base of pair to find (ex: 'G')\n" +
		"\toutfileName - name of final output file, which will contain region information from BED file and its counts of gained, lost, and constant CpG sites.\n\n" +
		"\tOPTIONS: \n" +
		"\t-compare\n" +
		"\tIn --compare mode, gained, lost, and constant pairs in genome 1 relative to genome 2 are counted.\n" +
		"\tAll paired base changes from substitutions, insertions, and deletions are reported.\n" +
		"\tWeak to strong (A/T -> C/G) and strong to weak (C/G -> A/T) substitutions are also reported.\n" +
		"\tNote: the algorithm counts seq1: C-G vs. seq2: CCG as a CG gain, followed by a CG loss, in seq1 relative to seq2.\n" +
		"\tSome may interpret this as constant CG site instead; check alignments if this is a concern.\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var compare *bool = flag.Bool("compare", false, "Will find paired base changes (gain, loss, cons) in sequence 1 relative to sequence 2 in the provided multiFas.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	//if number of provided args is not 3
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	//args
	s := Settings{
		FaDir:   flag.Arg(0),
		Bed:     flag.Arg(1),
		BaseOne: flag.Arg(2),
		BaseTwo: flag.Arg(3),
		Outfile: flag.Arg(4),
		Compare: *compare,
	}

	countByChrom(s)
}
