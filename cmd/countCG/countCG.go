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
)

type Settings struct {
	FaDir   string //input fasta directory with all chromosomes
	Bed     string //input bed file
	Outfile string //output file name
	Compare bool
}

type CpGcounts struct {
	Chrom string
	Name  string
	Gain  int
	Loss  int
	Cons  int
}

func isCons(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if firstSeq1 == dna.C && firstSeq2 == dna.G && secondSeq1 == dna.C && secondSeq2 == dna.G {
		return true
	} else {
		return false
	}
}

func isGain(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if firstSeq1 == dna.C && firstSeq2 == dna.G && !(secondSeq1 == dna.C && secondSeq2 == dna.G) {
		return true
	} else {
		return false
	}
}

func isLoss(firstSeq1, firstSeq2, secondSeq1, secondSeq2 dna.Base) bool {
	if !(firstSeq1 == dna.C && firstSeq2 == dna.G) && secondSeq1 == dna.C && secondSeq2 == dna.G {
		return true
	} else {
		return false
	}
}

func nextBase(region []dna.Base, currPos int) (nextBase dna.Base) {
	for i := currPos; i < len(region); i++ {
		if dna.DefineBase(region[i]) {
			return region[i]
		}
	}
	return dna.Gap
}

func findChromInFile(filedir string, chrom string) (string, error) {
	//09-29-25: force naming convention: chr1.fa, chr2.fa, etc. so chr name can be extracted easily
	files, err := os.ReadDir(filedir)
	if err != nil {
		return "", err
	}
	expectedName := chrom + ".fa"
	for _, f := range files {
		fileName := f.Name()
		if fileName == expectedName {
			filePath := filepath.Join(filedir, fileName)
			return filePath, nil
		}
	}
	return "", fmt.Errorf("File %s not found in %s", chrom+".fa", filedir)
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

	//add this region (chr, alignment start, alignment end, name) to slice of structs bedAln
	//newString := chrom + "\t" + strconv.Itoa(startAlnPos) + "\t" + strconv.Itoa(endAlnPos) + "\t" + name + "\n"
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

		//use the last reference/alignment position pair to convert the current region's ref start coordinate to an
		//alignment start coordinate
		startAlnPos = fasta.RefPosToAlnPosCounter(refSeq, region.ChromStart, lastRefPos, lastAlnPos)

		//use the same strategy for the current region's end coordinate
		endAlnPos = fasta.RefPosToAlnPosCounter(refSeq, region.ChromEnd, lastRefPos, lastAlnPos)

		//add this information to bedAln
		//newString := chrom + "\t" + strconv.Itoa(startAlnPos) + "\t" + strconv.Itoa(endAlnPos) + "\t" + name + "\n"
		output = append(output, bed.Bed{Chrom: chrom, ChromStart: startAlnPos, ChromEnd: endAlnPos, Name: regName, FieldsInitialized: 4})

		//reset the ref/alignment position pair to that of the current region, so that it can be used for the next region
		lastRefPos = region.ChromEnd
		lastAlnPos = endAlnPos
	}
	return output
}

func countCpG(fastaFile string, alnBed []bed.Bed) (regionCpGs []bed.Bed, err error) {

	chr := fasta.Read(fastaFile)

	//check that there are not more than 1 fasta sequence
	if len(chr) != 1 {
		log.Fatalf("Error: expecting exactly one record in fasta file, but got %d. If you want to compare 2 sequences, use --compare mode.\n", len(chr))
	}

	seq := chr[0].Seq
	//if sequence is empty
	if len(seq) == 0 {
		log.Fatalf("Error: fasta sequence is empty.\n")
	}

	for _, region := range alnBed {
		var f1, f2 dna.Base
		//get the following information
		chrom := region.Chrom
		start := region.ChromStart
		end := region.ChromEnd
		name := region.Name

		//Find the current region in the fasta sequences:
		//(Here, we assume that both the genome of interest and the ancestral genome are unmodified
		//from the original multiple alignment, so the same alignment coordinates apply to both)
		regSeq := seq[start:end]

		//initialize number of each category of CpG sites to 0

		cpgCount := 0

		//for length of entire sequence defined by alignment start/end coordinates
		for i := 0; i < len(regSeq)-1; i++ {
			f1, f2 = regSeq[i], regSeq[i+1]
			if f1 == dna.C && f2 == dna.G {
				cpgCount++
			}
		}
		regionCpGs = append(regionCpGs, bed.Bed{Chrom: chrom, ChromStart: start, ChromEnd: end, Name: name, Score: cpgCount, FieldsInitialized: 5})
	}
	return regionCpGs, nil
}

func compareCpGcount(fastaFile string, alnBed []bed.Bed) (regionCpGs []CpGcounts, err error) {
	seqs := fasta.Read(fastaFile)

	if len(seqs) != 2 {
		log.Fatalf("Error: expecting exactly two records in fasta file, but got %d. --compare mode compares exactly 2 aligned sequences.\n", len(seqs))
	}
	//define genome of interest and ancestral genome in the multiFa
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
		//(Here, we assume that both the genome of interest and the ancestral genome are unmodified
		//from the original multiple alignment, so the same alignment coordinates apply to both)
		firstRegseq := firstSeq[start:end]
		secondRegseq := secondSeq[start:end]

		//initialize number of each category of CpG sites to 0

		gainCount := 0
		lossCount := 0
		consCount := 0

		//for length of entire sequence defined by alignment start/end coordinates
		for i := 0; i < len(firstRegseq)-1; i++ {
			//f1, f2 and s1, s2 are calls to nextBase()
			f1, s1 := firstRegseq[i], secondRegseq[i]
			f2 := nextBase(firstRegseq, i+1)
			s2 := nextBase(secondRegseq, i+1)

			switch {
			case isCons(f1, f2, s1, s2):
				consCount++
			case isGain(f1, f2, s1, s2):
				gainCount++
			case isLoss(f1, f2, s1, s2):
				lossCount++
			}
		}
		regionCpGs = append(regionCpGs, CpGcounts{Chrom: chrom, Name: name, Gain: gainCount, Loss: lossCount, Cons: consCount})
	}
	return regionCpGs, nil
}

func countByChrom(s Settings) {
	var counts bed.Bed
	var countsComp CpGcounts

	bedBychrom := make(map[string][]bed.Bed)
	bedFile := s.Bed
	//preserve the original order of bed regions in the file
	bedRegions := bed.Read(bedFile)

	//make map by chromosome, to process one chromosome at a time (only open one fasta file at a time, loading it only once)
	for _, region := range bedRegions {
		bedBychrom[region.Chrom] = append(bedBychrom[region.Chrom], region)
	}

	file := fileio.EasyCreate(s.Outfile)

	if !s.Compare {
		fileio.WriteToFileHandle(file, "Chrom\tStart\tEnd\tName\tCpGcount")
	} else {
		fileio.WriteToFileHandle(file, "Chrom\tStart\tEnd\tName\tGain\tLoss\tCons")
	}

	//fast lookup table, PER CHROMOSOME, to match CpGcounts to the region's name in comp mode
	chromCountsMap := make(map[string]map[string]bed.Bed)
	//fast lookup table, PER CHROMOSOME, to match CpGcounts to the region's name in comp mode
	chromCountsCompMap := make(map[string]map[string]CpGcounts)

	for chrom, regions := range bedBychrom {
		fastaFile, err := findChromInFile(s.FaDir, chrom)
		exception.PanicOnErr(err)

		alnBed := RefPosToAlnPosBed(regions, fastaFile)

		if !s.Compare {
			counts, err := countCpG(fastaFile, alnBed)
			exception.PanicOnErr(err)

			cpgMap := make(map[string]bed.Bed)

			for _, r := range counts {
				cpgMap[r.Name] = r
			}
			chromCountsMap[chrom] = cpgMap
		} else {
			counts, err := compareCpGcount(fastaFile, alnBed)
			exception.PanicOnErr(err)

			cpgMap := make(map[string]CpGcounts)

			for _, c := range counts {
				cpgMap[c.Name] = c
			}
			chromCountsCompMap[chrom] = cpgMap
		}
	}

	// output regions in original BED order by looping through the regions from the original BED file
	for _, region := range bedRegions {
		if !s.Compare {
			countsMap, ok := chromCountsMap[region.Chrom]
			if !ok {
				log.Fatalf("No counts found for chromosome %s", region.Chrom)
			}
			counts, ok = countsMap[region.Name]
			if !ok {
				log.Fatalf("Counts not found for region %s", region.Name)
			}
		} else {
			countsCompMap, ok := chromCountsCompMap[region.Chrom]
			if !ok {
				log.Fatalf("No counts found for chromosome %s", region.Chrom)
			}
			countsComp, ok = countsCompMap[region.Name]
			if !ok {
				log.Fatalf("Counts not found for region %s", region.Name)
			}
		}

		if !s.Compare {
			line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d",
				region.Chrom, region.ChromStart, region.ChromEnd, region.Name, counts.Score)
			fileio.WriteToFileHandle(file, line)
		} else {
			line := fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%d",
				region.Chrom, region.ChromStart, region.ChromEnd, region.Name,
				countsComp.Gain, countsComp.Loss, countsComp.Cons)
			fileio.WriteToFileHandle(file, line)
		}
	}

	// Close file
	err := file.Close()
	exception.PanicOnErr(err)

	fmt.Println("CG counts found and written to", s.Outfile)
}

func usage() {
	fmt.Print("countCG fastaDir bedfileName outfileName\n" +
		"\tCounts gained, lost, and constant CpG sites within BED regions in genome 1 relative to (aligned) genome 2.\n" +
		"\tAll CpG changes from substitutions, insertions, and deletions are reported.\n\n" +
		"\tARGS: \n" +
		"\tfastaDir - path to directory containing multiFas for all chromosomes to analyze. Files must be named: {chromosome name}.fa (files can be gzipped).\n" +
		"\tbedfileName - name of BED file containing all genomic regions in which to count CpG changes. All regions are required to have a name in column 4.\n" +
		"\t**Note: please ensure that all chromosomes in the BED file have corresponding multiFas in the provided fasta directory.**\n" +
		"\toutfileName - name of final output file, which will contain region information from BED file and its counts of gained, lost, and constant CpG sites.\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var compare *bool = flag.Bool("compare", false, "Will find CpG changes (gain, loss, cons) in genome 1 relative to genome 2 in the provided multiFas.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	//if number of provided args is not 3
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		FaDir:   flag.Arg(0),
		Bed:     flag.Arg(1),
		Outfile: flag.Arg(2),
		Compare: *compare,
	}

	countByChrom(s)
}
