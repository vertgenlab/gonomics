package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
	"sort"
)

// writeMap writes out a pseudobulk map to an io.writer
func writeMap(mp map[string]float64, writer *fileio.EasyWriter) {
	var total float64
	var write string
	var writeSlice []string
	for i := range mp {
		total, _ = mp[i]
		write = fmt.Sprintf("%s\t%f", i, total)
		writeSlice = append(writeSlice, write)
	}
	sort.Strings(writeSlice)
	fileio.WriteToFileHandle(writer, "construct\tcounts")
	for _, i := range writeSlice {
		fileio.WriteToFileHandle(writer, i)
	}
}

// parseBam takes in a cellranger count bam file and pulls out reads that are representative of the UMI and also returns the construct associated with the UMI.
func parseBam(s starrSeq.ScStarrSeqSettings) {
	var k int = 0
	var constructName, umiBx, cellType string
	var count float64
	var found bool
	var bit int32
	var sc bool = false
	var populateUmiStruct bool = false
	var allConstructs, umiBxSlice, allCellTypes []string
	var umi starrSeq.UMI
	var umiSlice []starrSeq.UMI
	var err error
	var out, outSam *fileio.EasyWriter

	cellTypeMap := make(map[string]string)
	if s.ScAnalysis != "" { //create a bool for cellTypeAnalysis
		sc = true
		ck := starrSeq.ReadClusterKey(s.ScAnalysis)
		for _, i := range ck {
			cellTypeMap[i.Bx] = i.Cluster
			allCellTypes = append(allCellTypes, i.Cluster)
		}
	}

	if s.ByCell != "" || sc || s.BinCells > 0 || s.DetermineBins != "" || s.CountMatrix != "" {
		populateUmiStruct = true
	}

	ch, head := sam.GoReadToChan(s.InFile)

	if !s.NoOut {
		out = fileio.EasyCreate(s.OutFile)
	}

	if s.SamOut != "" {
		outSam = fileio.EasyCreate(s.SamOut)
		sam.WriteHeaderToFileHandle(outSam, head)
	}

	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	if s.Bed != "" {
		bedEntries := bed.Read(s.Bed)
		for _, i := range bedEntries {
			tree = append(tree, i)
		}
		selectTree = interval.BuildTree(tree)
	}

	pseudobulkMap := make(map[string]float64)

	for i := range ch { //iterate of over the chanel of sam.Sam
		if s.UmiSat != "" {
			readUmi, _, _ := sam.QueryTag(i, "UB")
			cb, _, _ := sam.QueryTag(i, "CB")
			umiBx = fmt.Sprintf("%s_%s", readUmi, cb)
			umiBxSlice = append(umiBxSlice, umiBx)
		}
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
		bit = num.(int32)
		if bit&8 == 8 { // bit 8 is the flag for a UMI that was used in final count. I call these "valid" UMIs.
			k++
			if s.Bed != "" {
				overlappedBedEntry := interval.Query(selectTree, i, "any")
				if len(overlappedBedEntry) != 1 {
					continue
				}
				constructName = overlappedBedEntry[0].(bed.Bed).Name
			} else {
				construct, _, _ := sam.QueryTag(i, "GX") // get the construct associated with valid UMI
				constructName = construct.(string)
			}
			if populateUmiStruct {
				//rand.Seed(s.setSeed)
				cell, _, _ := sam.QueryTag(i, "CB") // get the cell barcode associated with valid UMI
				cellType, found = cellTypeMap[cell.(string)]
				if found {
					umi = starrSeq.UMI{Bx: cell.(string), Cluster: cellType, Construct: constructName}
				} else {
					umi = starrSeq.UMI{Bx: cell.(string), Cluster: "undefined", Construct: constructName}
				}
				allConstructs = append(allConstructs, umi.Construct) //list of all constructs found in the bam
				umiSlice = append(umiSlice, umi)
			}
			if s.SamOut != "" { //write output as sam
				sam.WriteToFileHandle(outSam, i)
			}
			//pseudobulk, default behavior. Use a map [constructName]rawCount to add up reads per construct
			count, found = pseudobulkMap[constructName]
			if !found {
				pseudobulkMap[constructName] = 1
			} else {
				pseudobulkMap[constructName] = count + 1
			}
		}
	}

	fmt.Println("Found this many valid UMIs: ", k)

	if s.ByCell != "" {
		outByCell := fileio.EasyCreate(s.ByCell)
		for _, i := range umiSlice {
			fileio.WriteToFileHandle(outByCell, fmt.Sprintf("%s\t%s", i.Bx, i.Construct))
		}
		err = outByCell.Close()
		exception.PanicOnErr(err)
	}
	if s.UmiSat != "" {
		starrSeq.UmiSaturation(umiBxSlice, s.UmiSat)
	}
	if s.CountMatrix != "" {
		starrSeq.MakeCountMatrix(s, umiSlice, allConstructs)
	}
	if sc && !s.NoOut { //singleCellAnalysis, handles all the writing as well
		starrSeq.SingleCellAnalysis(s, umiSlice, allConstructs, allCellTypes, cellTypeMap, out)
	}
	if s.DetermineBins != "" {
		s.BinCells = starrSeq.DetermineIdealBins(s, umiSlice)
		fmt.Printf("Using %d bins\n", s.BinCells)
	}
	if s.BinCells > 0 { //goes to binning and pseudobulk. Handles normalization and writing as well
		binnedCells := starrSeq.DistributeCells(s.BinCells, umiSlice, false)
		starrSeq.BinnedPseudobulk(binnedCells, out, s.InputNormalize)
	}
	if s.InputNormalize != "" && !s.NoOut { //normalize pseudobulk
		starrSeq.InputNormalize(pseudobulkMap, s.InputNormalize)
	}
	if !s.NoOut && s.ScAnalysis == "" && s.BinCells < 1 { //write out pseudobulk
		writeMap(pseudobulkMap, out)
	}
	if !s.NoOut {
		err = out.Close()
		exception.PanicOnErr(err)
	}
	if s.SamOut != "" {
		err = outSam.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print("cellrangerBam -- Takes in a cellranger bam file of STARR-seq reads and parses the extra flags field to pull out the " +
		"representative read for each UMI and which construct it belongs to. Multiple GEM wells from the same STARR-seq experiment can be provided in a comma-delimited list " +
		"in the 'inFile' field. The output is a tab-delimited table of read-counts for each constructs.\n" +
		"NOTE: The default behavior of this function works best with STARR-seq libraries where constructs don't have much similarity with each other.\n" +
		"For libraries that need barcoding (like GWAS or cross-species comparisons) use the -bed option with a bed file corresponding to barcode regions.\n" +
		"Usage: \n" +
		"cellrangerBam [options] inFile outFile\n\n")
	flag.PrintDefaults()
}

func main() {
	var byCell *string = flag.String("byCell", "", "Report the construct that each UMI belongs to and which cell in which it was found in a tab-delimited table."+
		" The table will be sent to the provided file name.")
	var inputNorm *string = flag.String("inputNorm", "", "Takes in a tab delimited table with construct name and input normalization value")
	var samOut *string = flag.String("samOut", "", "Filter the input bam file by reads that have valid UMIs and send the output, in sam format, to the provdided file name")
	var cellTypeAnalysis *string = flag.String("cellTypeAnalysis", "", "Takes in a tab delimited file that has cell barcode and cell type identification. "+
		"The ouptut of options will be a matrix that has counts for each construct in each cell type. The Seurat command WhichCells() can be used to generate the required list.")
	var binCells *int = flag.Int("binCells", 0, "Number of bins to randomly assign cells to. The output will be a psudobulk table for each bin")
	var umiSat *string = flag.String("umiSat", "", "Create a UMI saturation curve of all reads in the input bam file. The table containing total reads count and UMI count"+
		" will be send to the provided file name")
	var transfectionNorm *string = flag.String("transfectionNorm", "", "Bam file from the same scSTARR-seq experiment containing cellranger count alignments to a transfection reporter reference"+
		" for cell cluster GFP normalization. Multiple bam files from different GEM wells can be provided in a comma-separated list")
	var bed *string = flag.String("bed", "", "Use a bed file for assigning reads to constructs instead of the GTF that was used in cellranger mkref. The bed file must have the "+
		"name of the construct in the fourth field. Recommended for constructs with barcodes.")
	var ncNorm *string = flag.String("ncNorm", "", "Reports single cell data as the ratio of reads for each construct to negative control reads in the same cell type."+
		"A file that contains a line-delimited list of negative control construct names must be provided.")
	var determineBins *string = flag.String("determineBins", "", "Determine the ideal number of pseudobulk bins to use. "+
		"The ideal number of bins will be determined to be the number where all bins have at least 1 read from all negative control constructs. "+
		"A file that contains a line-delimited list of negative control construct names must be provided.")
	var inputSequencing *string = flag.String("inputSequencing", "", "Returns library statistics (Reads per 500bp, percent of each construct in library, and input normalization factor"+
		"for an input library sequencing sam/bam file. A bed file corresponding to constructs on the reference that the input file is aligned to must be provided and the name of the construct must"+
		"be in the 4th field. The input alignment file must be sorted by position. PCR duplicates will be removed.")
	//var setSeed *int64 = flag.Int64("setSeed", 1, "Set a seed for random number generation.")
	var countMatrix *string = flag.String("scCount", "", "Create count matrix that has cell barcode in rows and constructs as columns. If GFP bam files are provided with the"+
		" -gfpNorm option, GFP will be included as a column. The count matrix will be sent to the file name provided")
	var noOut *bool = flag.Bool("noOut", false, "If the -noOut option is used a psuedobulk table will not be created a the function will only take in the inFile argument.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	var expectedNumArgs int

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
	if *binCells > 0 && *cellTypeAnalysis != "" {
		log.Fatalf("Bin cells cannot be used with -cellTypeAnalysis")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	var s starrSeq.ScStarrSeqSettings = starrSeq.ScStarrSeqSettings{
		InFile:           flag.Arg(0),
		OutFile:          flag.Arg(1),
		InputNormalize:   *inputNorm,
		ByCell:           *byCell,
		ScAnalysis:       *cellTypeAnalysis,
		BinCells:         *binCells,
		UmiSat:           *umiSat,
		SamOut:           *samOut,
		TransfectionNorm: *transfectionNorm,
		Bed:              *bed,
		NcNorm:           *ncNorm,
		DetermineBins:    *determineBins,
		InputSequencing:  *inputSequencing,
		//setSeed:         *setSeed,
		CountMatrix: *countMatrix,
		NoOut:       *noOut,
	}

	//deal with multiple gem wells
	var multipleGems bool
	var multipleGfpGems bool = false
	var err error

	s.InFile, multipleGems = starrSeq.DetectMultipleGems(s.InFile, s.ScAnalysis)
	if *transfectionNorm != "" {
		s.TransfectionNorm, multipleGfpGems = starrSeq.DetectMultipleGfpGems(s.TransfectionNorm)
	}

	//output or input sequencing
	if *inputSequencing != "" {
		starrSeq.ParseInputSequencingSam(s)
	} else {
		parseBam(s)
	}

	//if multiple gem wells were used remove the tmp files
	if multipleGems {
		err = os.Remove(s.InFile)
		exception.PanicOnErr(err)
	}
	if multipleGfpGems {
		err = os.Remove(s.TransfectionNorm)
		exception.PanicOnErr(err)
	}
}
