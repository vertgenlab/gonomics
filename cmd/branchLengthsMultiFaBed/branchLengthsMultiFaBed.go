package main

import (
	"fmt"
	"flag"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/phylo"
	"log"
)

type Settings struct {
	Chrom string
	InFaFile string
	InBedFile string
	VelLengthBedFile string
	InitialLengthBedFile string
	NumUngappedSitesBedFile string
	SearchSpaceBed string
	SearchSpaceProportion float64
	UseSnpDistance bool
	Verbose bool
	Epsilon float64
	AllowNegative bool
	ZeroDistanceWeightConstant float64
}

func branchLengthsMultiFaBed(s Settings) {
	records := fasta.Read(s.InFaFile)
	var bitArray []bool
	var reachedEnd bool = false
	var currLengths phylo.AccelBranchLengths
	var currSize, currAln, currRef, currUngappedCount int = 0, 0, 0, 0
	var currDistance phylo.AccelDistances = phylo.AccelDistances{0, 0, 0, 0, 0, 0}
	var err error

	//a couple of safety checks
	if len(records) != 4 {
		log.Fatalf("brancLneghtsMultiFaBed accepts a multiFa file with 4 records, found %v.", len(records))
	}
	if len(records[1].Seq) != len(records[0].Seq) || len(records[2].Seq) != len(records[0].Seq) || len(records[3].Seq) != len(records[0].Seq) {
		log.Fatalf("Error. All records must be of the same sequence length.")
	}

	referenceLength := fasta.AlnPosToRefPos(records[0], len(records[0].Seq) - 1)
	if s.SearchSpaceBed != "" {
		bitArray = phylo.MakeBitArrayFromSearchSpaceBed(s.SearchSpaceBed, referenceLength, s.Chrom)
	}

	regions := bed.Read(s.InBedFile)
	bed.SortByCoord(regions)
	VelOut := fileio.EasyCreate(s.VelLengthBedFile)
	InitialOut := fileio.EasyCreate(s.InitialLengthBedFile)
	UngappedBasesOut := fileio.EasyCreate(s.NumUngappedSitesBedFile)

	for i := 0; i < len(regions); i++ {
		if passesSearchSpaceTest(regions[i], bitArray, s) {
			currSize = regions[i].ChromEnd - regions[i].ChromStart
			currAln = fasta.RefPosToAlnPosCounter(records[0], regions[i].ChromStart, currRef, currAln)//we update the currAln to the position of the next bed entry, using the cached ref and aln from the last iteration
			currRef = regions[i].ChromStart//now we can safely update currRef to the position of the current bed entry.
			if s.UseSnpDistance {
				reachedEnd = phylo.AccelFourWaySnpDistances(records, currAln, currSize, &currDistance)
			} else {
				reachedEnd = phylo.AccelFourWayMutationDistances(records, currAln, currSize, &currDistance)
			}
			if reachedEnd {
				log.Fatalf("Error: bed entry ran off the end of the multiple alignment chromosome. Offending entry starts at position %s\t%v", regions[i].Chrom, regions[i].ChromStart)
			}
			currLengths = phylo.BranchLengthsAlternatingLeastSquares(currDistance, s.AllowNegative, s.Verbose, s.ZeroDistanceWeightConstant, s.Epsilon)
			currUngappedCount = numUngappedInBedRange(records, currAln, currSize)

			bed.WriteBed(VelOut, bed.Bed{Chrom: s.Chrom, ChromStart: regions[i].ChromStart, ChromEnd: regions[i].ChromEnd, Name: fmt.Sprintf("%v", currLengths.B1), FieldsInitialized: 4})
			bed.WriteBed(InitialOut, bed.Bed{Chrom: s.Chrom, ChromStart: regions[i].ChromStart, ChromEnd: regions[i].ChromEnd, Name: fmt.Sprintf("%v", currLengths.B3), FieldsInitialized: 4})
			bed.WriteBed(UngappedBasesOut, bed.Bed{Chrom: s.Chrom, ChromStart: regions[i].ChromStart, ChromEnd: regions[i].ChromEnd, Name: fmt.Sprintf("%v", currUngappedCount), FieldsInitialized: 4})
		}
	}

	err = VelOut.Close()
	exception.PanicOnErr(err)
	err = InitialOut.Close()
	exception.PanicOnErr(err)
	err = UngappedBasesOut.Close()
	exception.PanicOnErr(err)
}

func numUngappedInBedRange(records []fasta.Fasta, currAln int, currSize int) int {
	var baseCount, i, unGappedCount int = 0, 0, 0
	for i = currAln; baseCount < currSize && i < len(records[0].Seq); i++ {
		if phylo.IsUngappedColumn(records, i) {
			unGappedCount++
		}
	}
	return unGappedCount
}

func passesSearchSpaceTest(b bed.Bed, bitArray []bool, s Settings) bool {
	if s.SearchSpaceBed == "" {
		return true
	}
	if b.Chrom != s.Chrom {
		return false
	}
	var count int = 0
	for i := b.ChromStart; i < b.ChromEnd; i++ {
		if bitArray[i] {
			count++
		}
	}
	proportion := float64(count) / float64(b.ChromEnd - b.ChromStart)
	return proportion >= s.SearchSpaceProportion
}

func usage() {
	fmt.Print(
		"branchLengthsMultiFaBed - Using a four-way multiple alignment (including the reference species followed by three successive outgroups),\n" +
			"this program calculates branch lengths in units of estimate substitutions for regions in the alignment specified by an input bed file.\n" +
			"Two branches are returned: the branch for the 'vel' branch which separates the reference from\n" +
			"the common ancestor with the first outspecies, and the 'initial' branch length, which separates the first common ancestor from the common ancestor with the second\n" +
			"outspecies. Uses the Fitch-Margoliash method for branch length estimation. \n" +
			"Usage:\n" +
			"branchLengthsMultiFaBed chromName in.fa in.bed velLength.bed initialLength.bed NumUngappedSites.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 6
	var searchSpaceBed *string = flag.String("searchSpaceBed", "", "This bed file should contain the valid alignment regions within the multiFa. Input bed entries that do not overlap the search space will not appear in the output file.")
	var searchSpaceProportion *float64 = flag.Float64("searchSpaceProportion", 0.5, "Proportion of bed entry that must overlap search space in order to be evaluated.")
	var useSnpDistance *bool = flag.Bool("useSnpDistance", false, "Calculate pairwise distances with SNPs instead of the default mutation distance, which counts INDELs.")
	var verbose *bool = flag.Bool("verbose", false, "Enables debug prints.")
	var epsilon *float64 = flag.Float64("epsilon", 1e-8, "Set the error threshold for alternating least squares branch length calculation.")
	var allowNegative *bool = flag.Bool("allowNegative", false, "Allow the algorithm to evaluate negative branch lengths. This program will constrain the optimal solution to non-negative branch lengths by default.")
	var zeroDistanceWeightConstant *float64 = flag.Float64("zeroDistanceWeightConstant", 1000, "Set the relative error weight applied to pairs of species with a pairwise distance of zero.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	chromName := flag.Arg(0)
	inFaFile := flag.Arg(1)
	inBedFile := flag.Arg(2)
	veLengthBedFile := flag.Arg(3)
	initialLengthBedFile := flag.Arg(4)
	numUngappedSitesBedFile := flag.Arg(5)

	s := Settings {
		Chrom: chromName,
		InFaFile: inFaFile,
		InBedFile: inBedFile,
		VelLengthBedFile: veLengthBedFile,
		InitialLengthBedFile: initialLengthBedFile,
		NumUngappedSitesBedFile: numUngappedSitesBedFile,
		SearchSpaceBed: *searchSpaceBed,
		SearchSpaceProportion: *searchSpaceProportion,
		UseSnpDistance: *useSnpDistance,
		Verbose: *verbose,
		Epsilon: *epsilon,
		AllowNegative: *allowNegative,
		ZeroDistanceWeightConstant: *zeroDistanceWeightConstant,
	}

	branchLengthsMultiFaBed(s)
}