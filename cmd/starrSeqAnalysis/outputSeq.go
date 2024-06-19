package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
)

func outputSeqUsage() {

}

func parseOutputSeqArgs() {
	var expectedNumArgs int

	outputSeqFlags := flag.NewFlagSet("outputSeq", flag.ExitOnError)

	var validUmis *string = outputSeqFlags.String("validUmis", "", "Report the cell barcode, UMI, construct and cell type, if applicable, for each valid / unique UMI."+
		" The table will be sent to the provided file name.")
	var inputNorm *string = outputSeqFlags.String("inputNorm", "", "Takes in a tab delimited table with construct name and input normalization value")
	var samOut *string = outputSeqFlags.String("samOut", "", "Filter the input bam file by reads that have valid UMIs and send the output, in sam format, to the provdided file name")
	var cellTypeAnalysis *string = outputSeqFlags.String("cellTypeAnalysis", "", "Takes in a tab delimited file that has cell barcode and cell type identification. "+
		"The ouptut of options will be a matrix that has counts for each construct in each cell type. The Seurat command WhichCells() can be used to generate the required list.")
	var binCells *int = outputSeqFlags.Int("binCells", 0, "Number of bins to randomly assign cells to. The output will be a psudobulk table for each bin")
	var umiSat *string = outputSeqFlags.String("umiSat", "", "Create a UMI saturation curve of all reads in the input bam file. The table containing total reads count and UMI count"+
		" will be send to the provided file name")
	var transfectionNorm *string = outputSeqFlags.String("transfectionNorm", "", "Bam file from the same scSTARR-seq experiment containing cellranger count alignments to a transfection reporter reference"+
		" for cell cluster GFP normalization. Multiple bam files from different GEM wells can be provided in a comma-separated list")
	var bed *string = outputSeqFlags.String("bed", "", "Use a bed file for assigning reads to constructs instead of the GTF that was used in cellranger mkref. The bed file must have the "+
		"name of the construct in the fourth field. Recommended for constructs with barcodes.")
	var ncNorm *string = outputSeqFlags.String("ncNorm", "", "Reports single cell data as the ratio of reads for each construct to negative control reads in the same cell type."+
		"A file that contains a line-delimited list of negative control construct names must be provided.")
	var determineBins *string = outputSeqFlags.String("determineBins", "", "Determine the ideal number of pseudobulk bins to use. "+
		"The ideal number of bins will be determined to be the number where all bins have at least 1 read from all negative control constructs. "+
		"A file that contains a line-delimited list of negative control construct names must be provided.")
	//var setSeed *int64 = flag.Int64("setSeed", 1, "Set a seed for random number generation.")
	var countMatrix *string = outputSeqFlags.String("countMatrix", "", "Create count matrix that has cell barcode in rows and constructs as columns. If GFP bam files are provided with the"+
		" -transfectionNorm option, GFP will be included as a column. If -countMatrixCellTypes is used, an additional column for cell type will be included. "+
		"The count matrix will be sent to the file name provided")
	var noOut *bool = outputSeqFlags.Bool("noOut", false, "If the -noOut option is used a psuedobulk table will not be created a the function will only take in the inFile argument.")
	var altMapping *string = outputSeqFlags.String("altMapping", "", "Use a bed file for construct mapping and collapse corrected-UMI by hand rather than rely on cellranger count valid UMIs. "+
		"Provde a bed file with construct names for construct determination.")
	var countMatrixCellTypes *string = outputSeqFlags.String("countMatrixCellTypes", "", "Provide a tab-delimited file of cellBarcode -- cell type to add an additional column"+
		"to a count matrix corresponding to cell type.")
	var stats *string = flag.String("stats", "", "Provide a line-delimited list of negative controll constructs that all constructs will be compared against with a "+
		"Wilcoxin rank sum test. An additional column with pValue will be added to the pseudobulk map. Must be used with -binCells or -determineBins.")

	err := outputSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	outputSeqFlags.Usage = func() { makeRefUsage(outputSeqFlags) }

	if *noOut {
		expectedNumArgs = 1
	} else {
		expectedNumArgs = 2
	}

	if *ncNorm != "" && *transfectionNorm != "" {
		log.Fatalf("Error: only one normalization method for single-cell data can be used at once")
	}
	if *binCells < 0 {
		log.Fatalf("Error: -binCells must be a positive intiger")
	}
	if (*binCells > 0 || *determineBins != "") && *cellTypeAnalysis != "" {
		log.Fatalf("Bin cells cannot be used with -cellTypeAnalysis")
	}
	if *altMapping != "" && *samOut != "" {
		log.Fatalf("altMapping is not compatable with samOut. If you are interested in collapsing UMI in a sam file use samFilter -collapseUMI")
	}
	if *stats != "" && (*binCells < 2 || *determineBins != "") {
		log.Fatalf("-stats must be used with -binCells or -determineBins")
	}

	s := starrSeq.OutputSeqSettings{
		InFile:           flag.Arg(0),
		OutFile:          flag.Arg(1),
		InputNormalize:   *inputNorm,
		ValidUmis:        *validUmis,
		ScAnalysis:       *cellTypeAnalysis,
		BinCells:         *binCells,
		UmiSat:           *umiSat,
		SamOut:           *samOut,
		TransfectionNorm: *transfectionNorm,
		Bed:              *bed,
		NcNorm:           *ncNorm,
		DetermineBins:    *determineBins,
		//setSeed:         *setSeed,
		CountMatrix:          *countMatrix,
		NoOut:                *noOut,
		AltMapping:           *altMapping,
		CountMatrixCellTypes: *countMatrixCellTypes,
		Stats:                *stats,
	}

	if len(outputSeqFlags.Args()) != expectedNumArgs {
		outputSeqFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(outputSeqFlags.Args()))
	}
	if *altMapping != "" {
		starrSeq.Alt(s)
	} else {
		starrSeq.ParseCellrangerBam(s)
	}
}
