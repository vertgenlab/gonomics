package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type Settings struct {
	InFa    string //input fasta file (1 chromosome)
	Chrom   string //chrom to analyze
	BaseOne string //base one in pair
	BaseTwo string //base two in pair
	Outfile string //output file name
	BedFile string //input bed file
	Compare bool   //-compare mode for 2 aligned sequences
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

// RefPosToAlnPosBed loops existing functions for converting reference positions to alignment positions
// (RefPosToAlnPos and RefPostoAlnPosCount), to convert the coordinates of an entire bed file.
// It takes a bed file and a multiFa as inputs and outputs the bed regions with alignment coordinates.
func RefPosToAlnPosBed(input []bed.Bed, multiFa []fasta.Fasta) (output []bed.Bed) {
	var regName string

	refSeq := multiFa[0]
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
func countPairOfBasesHelper(inFa fasta.Fasta, b1, b2 dna.Base) (counts int) {
	var f1, f2 dna.Base

	seq := inFa.Seq

	//if sequence is empty
	if len(seq) == 0 {
		log.Fatalf("Error: fasta sequence is empty.\n")
	}

	//initialize number of paired sites to 0
	pairOfBaseCount := 0

	//for length of entire sequence defined by reference start/end coordinates
	for i := 0; i < len(seq)-1; i++ {
		f1, f2 = seq[i], seq[i+1]
		if f1 == b1 && f2 == b2 {
			pairOfBaseCount++
		}
	}
	//output
	return pairOfBaseCount
}

// -compare mode counting gains, losses, and constant pairs of bases provided by user
// inputs: multiFasta for the correct chromosome returned by findChromInFile, alignment-converted BED returned by RefPosToAlnPosBed, base 1 and base 2 to count as a pair
func comparePairOfBaseCount(firstFa, secondFa fasta.Fasta, b1, b2 dna.Base) (gainCount, lossCount, consCount, weakToStrongCount, strongToWeakCount int) {
	//define sequence 1 and sequence 2 in the multiFa
	firstSeq := firstFa.Seq
	secondSeq := secondFa.Seq

	//initialize number of each category of paired sites to 0
	gainCount = 0
	lossCount = 0
	consCount = 0

	//initialize each category of substitutions to 0
	weakToStrongCount = 0
	strongToWeakCount = 0

	//for length of entire sequence defined by alignment start/end coordinates
	for i := 0; i < len(firstSeq)-1; i++ {
		//f1, f2 and s1, s2 are calls to nextBase()
		f1, s1 := firstSeq[i], secondSeq[i]
		f2 := nextBase(firstSeq, i+1)
		s2 := nextBase(secondSeq, i+1)

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
		if i == len(firstSeq)-2 {
			if isWeakToStrong(f2, s2) {
				weakToStrongCount++
			}

			if isStrongToWeak(f2, s2) {
				strongToWeakCount++
			}
		}
	}
	return gainCount, lossCount, consCount, weakToStrongCount, strongToWeakCount
}

// perform above^ functions to count pairs of bases in each FASTA or BED region
func countPairOfBases(s Settings) {
	//convert user-provided bases to dna.Base type
	if len(s.BaseOne) != 1 || len(s.BaseTwo) != 1 {
		log.Fatalf("Error: Enter one DNA base for 'base one' and one DNA base for 'base two'.\n")
	}

	baseOne := dna.StringToBase(strings.TrimSpace(s.BaseOne))
	baseTwo := dna.StringToBase(strings.TrimSpace(s.BaseTwo))

	inFasta := fasta.Read(s.InFa)

	outfile := fileio.EasyCreate(s.Outfile)

	if !s.Compare {
		if s.BedFile == "" {
			if len(inFasta) != 1 {
				log.Fatalf("Error: expecting exactly one record in fasta file, but got %d. If you want to compare 2 sequences, use --compare mode.\n", len(inFasta))
			}
			inFaEntry := inFasta[0]
			pairCounts := countPairOfBasesHelper(inFaEntry, baseOne, baseTwo)

			fileio.WriteToFileHandle(outfile, "Chrom\tPairOfBasesCount")
			fileio.WriteToFileHandle(outfile, fmt.Sprintf("%s\t%d",
				s.Chrom, pairCounts))
		} else {
			fileio.WriteToFileHandle(outfile, "Chrom\tStart\tEnd\tName\tPairOfBasesCount")
			bedRegions := bed.Read(s.BedFile)
			for _, region := range bedRegions {
				if region.Chrom != s.Chrom {
					log.Fatalf("Error: Chromosome in BED region does not match.")
				}
				inSubFa := convert.SingleBedToFasta(region, inFasta)
				regPairCounts := countPairOfBasesHelper(inSubFa, baseOne, baseTwo)
				line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d", region.Chrom, region.ChromStart, region.ChromEnd, region.Name, regPairCounts)
				fileio.WriteToFileHandle(outfile, line)
			}
		}
	} else {
		if s.BedFile == "" {
			if len(inFasta) != 2 {
				log.Fatalf("Error: expecting exactly two records in fasta file, but got %d. --compare mode compares exactly 2 aligned sequences.\n", len(inFasta))
			}
			seqOne := inFasta[0]
			seqTwo := inFasta[1]

			gain, loss, cons, weakToStrong, strongToWeak := comparePairOfBaseCount(seqOne, seqTwo, baseOne, baseTwo)
			fileio.WriteToFileHandle(outfile, "Chrom\tGain\tLoss\tCons\tWeakToStrongSubs\tStrongToWeakSubs")
			fileio.WriteToFileHandle(outfile, fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%d", s.Chrom, gain, loss, cons, weakToStrong, strongToWeak))
		} else {
			fileio.WriteToFileHandle(outfile, "Chrom\tStart\tEnd\tName\tGain\tLoss\tCons\tWeakToStrongSubs\tStrongToWeakSubs")
			bedRegions := bed.Read(s.BedFile)
			//make map to key by name for easy ref coordinate writing to file
			bedMap := make(map[string]bed.Bed)
			for _, region := range bedRegions {
				if region.Chrom != s.Chrom {
					log.Fatalf("Error: Chromosome in BED region does not match.")
				}
				bedMap[region.Name] = region
			}

			if len(inFasta) != 2 {
				log.Fatalf("Error: expecting exactly two records in fasta file, but got %d. --compare mode compares exactly 2 aligned sequences.\n", len(inFasta))
			}
			alnBed := RefPosToAlnPosBed(bedRegions, inFasta)
			for _, alnRegion := range alnBed {

				firstRegSeq := fasta.Extract(inFasta[0], alnRegion.ChromStart, alnRegion.ChromEnd, alnRegion.Name)
				secondRegSeq := fasta.Extract(inFasta[1], alnRegion.ChromStart, alnRegion.ChromEnd, alnRegion.Name)

				gain, loss, cons, weakToStrong, strongToWeak := comparePairOfBaseCount(firstRegSeq, secondRegSeq, baseOne, baseTwo)
				ref := bedMap[alnRegion.Name]
				fileio.WriteToFileHandle(outfile, fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d", ref.Chrom, ref.ChromStart, ref.ChromEnd, ref.Name, gain, loss, cons, weakToStrong, strongToWeak))
			}
		}
	}
	//close file
	err := outfile.Close()
	exception.PanicOnErr(err)

	fmt.Println("Pair counts found and written to", s.Outfile)
}

func usage() {
	fmt.Print("countPairOfBases [options] fastaDir chromName baseOne baseTwo outfileName\n" +
		"\tCounts pairs of bases (ex: CG) within BED regions.\n" +
		"\tARGS: \n" +
		"\tfastaFile - path to fasta/multiFasta for the chromosome to analyze (files can be gzipped).\n" +
		"\tchromName - name of chromosome to analyze in format: 'chr1', 'chrX', etc.\n" +
		"\tbaseOne - first base of pair to find (ex: 'C')\n" +
		"\tbaseTwo - second base of pair to find (ex: 'G')\n" +
		"\toutfileName - name of final output file, which will contain region information from BED file and its counts of gained, lost, and constant CpG sites.\n\n" +
		"\tOPTIONS: \n" +
		"\tbedfile - name of BED file containing all genomic regions in which to count CpG changes. All regions are required to have a name in column 4.\n" +
		"\t**Note: please ensure that BED file is split by chromosome/only contains regions on the chromosome provided in the fastaFile/chromName.**\n" +
		"\t-compare\n" +
		"\tIn --compare mode, fastaFile is a 2-entry multiFasta, and gained, lost, and constant pairs in genome 1 relative to genome 2 are counted.\n" +
		"\tAll paired base changes from substitutions, insertions, and deletions are reported.\n" +
		"\tWeak to strong (A/T -> C/G) and strong to weak (C/G -> A/T) substitutions are also reported.\n" +
		"\tNote: the algorithm counts seq1: C-G vs. seq2: CCG as a CG gain, followed by a CG loss, in seq1 relative to seq2.\n" +
		"\tSome may interpret this as constant CG site instead; check alignments if this is a concern.\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var bedFile *string = flag.String("bedFile", "", "BED file path. Returns counts for all regions specified by the input BED file. Can be used in default or --compare mode.")
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
		InFa:    flag.Arg(0),
		Chrom:   flag.Arg(1),
		BaseOne: flag.Arg(2),
		BaseTwo: flag.Arg(3),
		Outfile: flag.Arg(4),
		BedFile: *bedFile,
		Compare: *compare,
	}

	countPairOfBases(s)
}
