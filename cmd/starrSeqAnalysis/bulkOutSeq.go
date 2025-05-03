package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/starrSeq"
	"golang.org/x/exp/slices"
	"log"
	"os"
	"strings"
)

func bulkOutSeqUsage(bulkOutSeqFlags *flag.FlagSet) {
	fmt.Print("starrSeqAnalysis bulkOutput [options] in.sam in.bed out.counts\n" +
		"Multiple replicates of the same STARR-seq experiment can be analysed at once by providing a comma-delimited list of alignment files instead of a single alignment file.\n")
	bulkOutSeqFlags.PrintDefaults()
}

func writeResults(outFile string, res []string, s starrSeq.BulkOutputSeqSettings) {
	out := fileio.EasyCreate(outFile)
	base := "construct\trawCounts"
	if s.InputNorm != "" {
		base += "\tnormCounts\tnormFactor"
	}
	if s.IntronStats {
		base += "\tintronCounts\tintronPercent"
	}
	if s.ZScore != "" {
		base += "\tenhZScore"
	}
	base += "\trep"
	fileio.WriteToFileHandle(out, base)
	slices.Sort(res)
	for i := range res {
		fileio.WriteToFileHandle(out, res[i])
	}
	exception.PanicOnErr(out.Close())
}

func parseBulkOutputArgs() {
	var expectedNumArgs int = 3
	bulkOutSeqFlags := flag.NewFlagSet("bulkOutput", flag.ExitOnError)
	var normFactor *string = bulkOutSeqFlags.String("normFactor", "", "Provide the output of starrSeqAnalysis inputSeq for input normalization")
	var dualBx *bool = bulkOutSeqFlags.Bool("dualBx", false, "The bed file provided reflects barcodes that are on each side of both construct. The bed file must be formatted with each barcode having the name: \"constructName_a\" or \"constructName_b\"")
	var pe *bool = bulkOutSeqFlags.Bool("pe", false, "The output sequencing file is from paired-end sequencing and name sorted with samtools sort -n")
	var singleBx *bool = bulkOutSeqFlags.Bool("singleBx", false, "Use if the bed file contains a single barcode region for each construct. Important to use because the program will"+
		"assume that the sam read will have to completely encompass the barcode region")
	var zscore *string = bulkOutSeqFlags.String("zscore", "", "provide a line-delimited text file with the names of the constructs that all other constructs will be normalized to. The output file will have an extra column containing the z-score of that construct")
	var checkBx *string = bulkOutSeqFlags.String("checkBx", "", "Provide the reference genome used for the alignment to check that the alignment-overlapped barcode doesn't have any mismatches. Must be used with either -singleBx or -dualBx")
	var intronStats *bool = bulkOutSeqFlags.Bool("intronStats", false, "Add percent reads with an intron to the output file")

	err := bulkOutSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)

	bulkOutSeqFlags.Usage = func() { bulkOutSeqUsage(bulkOutSeqFlags) }

	if len(bulkOutSeqFlags.Args()) != expectedNumArgs {
		bulkOutSeqFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(bulkOutSeqFlags.Args()))
	}

	if *dualBx && *singleBx {
		log.Fatalf("Error: -singleBx and -dualBx are incompatable with each other\n")
	}

	if *checkBx != "" && (!*dualBx && !*singleBx) {
		log.Fatalf("Error: -checkBx must be used with either -singleBx or -dualBx")
	}

	inSams := strings.Split(bulkOutSeqFlags.Arg(0), ",")

	var s starrSeq.BulkOutputSeqSettings
	var res []string
	for i := range inSams {
		s = starrSeq.BulkOutputSeqSettings{
			InSam:       inSams[i],
			InBed:       bulkOutSeqFlags.Arg(1),
			OutCounts:   bulkOutSeqFlags.Arg(2),
			InputNorm:   *normFactor,
			DualBx:      *dualBx,
			PairedEnd:   *pe,
			SingleBx:    *singleBx,
			ZScore:      *zscore,
			RepNum:      i + 1,
			CheckBx:     *checkBx,
			IntronStats: *intronStats,
		}
		res = append(res, starrSeq.BulkOutputSeq(s)...)
	}
	writeResults(bulkOutSeqFlags.Arg(2), res, s)
}
