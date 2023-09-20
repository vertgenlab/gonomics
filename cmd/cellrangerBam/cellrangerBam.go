package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"strings"
)

// parseBam takes in a cellranger count bam file and pulls out reads that are representative of the UMI and also returns the construct associated with the UMI.
func parseCellrangerBam(s starrSeq.ScStarrSeqSettings) {
	var k int = 0
	var constructName, umiBx, cellType, cellFormat string
	var found bool
	var bit uint8
	var multipleGems = false
	var gemNumber int = 1
	var allConstructs, umiBxSlice, allCellTypes []string
	var read starrSeq.Read
	var readSlice []starrSeq.Read
	var outSam *fileio.EasyWriter

	//detect multiple bam files
	files := strings.Split(s.InFile, ",")
	if len(files) > 1 {
		multipleGems = true
		if s.ScAnalysis != "" {
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -gfpNorm. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		}
	}
	cellTypeMap := make(map[string]string)
	if s.ScAnalysis != "" {
		ck := starrSeq.ReadClusterKey(s.ScAnalysis)
		for _, i := range ck {
			cellTypeMap[i.Bx] = i.Cluster
			allCellTypes = append(allCellTypes, i.Cluster)
		}
	} else if s.CountMatrixCellTypes != "" {
		ck := starrSeq.ReadClusterKey(s.CountMatrixCellTypes)
		for _, i := range ck {
			cellTypeMap[i.Bx] = i.Cluster
			allCellTypes = append(allCellTypes, i.Cluster)
		}
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

	for _, bam := range files {
		ch, head := sam.GoReadToChan(bam)
		if s.SamOut != "" {
			if gemNumber == 1 {
				outSam = fileio.EasyCreate(s.SamOut)
				sam.WriteHeaderToFileHandle(outSam, head)
			}
		}
		for i := range ch { //iterate of over the chanel of sam.Sam
			if s.UmiSat != "" {
				readUmi, foundUMI, _ := sam.QueryTag(i, "UB")
				cb, foundCb, _ := sam.QueryTag(i, "CB")

				if foundCb && foundUMI {
					umiBx = fmt.Sprintf("%s_%s", readUmi, cb)
					umiBxSlice = append(umiBxSlice, umiBx)
				}
			}
			num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
			bit = num.(uint8)
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
				//rand.Seed(s.setSeed)
				cell, _, _ := sam.QueryTag(i, "CB") // get the cell barcode associated with valid UMI
				if multipleGems {
					cellFormat = fmt.Sprintf("%s_%d", cell.(string), gemNumber)
				} else {
					cellFormat = cell.(string)
				}
				readUmi, _, _ := sam.QueryTag(i, "UB") // get the UMI associated with the valid UMI
				cellType, found = cellTypeMap[cellFormat]
				switch {
				case (s.ScAnalysis != "" || s.CountMatrixCellTypes != "") && found: //we want cell type data and its found
					read = starrSeq.Read{Bx: cellFormat, Cluster: cellType, Construct: constructName, UMI: readUmi.(string)}
				case (s.ScAnalysis != "" || s.CountMatrixCellTypes != "") && !found: //we want cell type data, but we didn't find the cell type in the map
					read = starrSeq.Read{Bx: cellFormat, Cluster: "undefined", Construct: constructName, UMI: readUmi.(string)}
				case s.ScAnalysis == "" && s.CountMatrixCellTypes == "": // we don't care about cell type analysis. cluster is nil
					read = starrSeq.Read{Bx: cellFormat, Construct: constructName, UMI: readUmi.(string)}
				}
				allConstructs = append(allConstructs, read.Construct) //list of all constructs found in the bam
				readSlice = append(readSlice, read)
				if s.SamOut != "" { //write output as sam
					sam.WriteToFileHandle(outSam, i)
				}
			}
		}
		gemNumber++
	}
	fmt.Println("Found this many valid UMIs: ", k)
	r := starrSeq.ReadSliceAnalysisSettings{ReadSlice: readSlice, FuncSettings: s, AllCellTypes: allCellTypes, AllConstructs: allConstructs, CellTypeMap: cellTypeMap, UmiBxSlice: umiBxSlice}
	starrSeq.ReadSliceAnalysis(r)
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
	var countMatrix *string = flag.String("countMatrix", "", "Create count matrix that has cell barcode in rows and constructs as columns. If GFP bam files are provided with the"+
		" -transfectionNorm option, GFP will be included as a column. If -countMatrixCellTypes is used, an additional column for cell type will be included. "+
		"The count matrix will be sent to the file name provided")
	var noOut *bool = flag.Bool("noOut", false, "If the -noOut option is used a psuedobulk table will not be created a the function will only take in the inFile argument.")
	var altMapping *string = flag.String("altMapping", "", "Use a bed file for construct mapping and collapse corrected-UMI by hand rather than rely on cellranger count valid UMIs. "+
		"Provde a bed file with construct names for construct determination.")
	var countMatrixCellTypes *string = flag.String("countMatrixCellTypes", "", "Provide a tab-delimited file of cellBarcode -- cell type to add an additional column"+
		"to a count matrix corresponding to cell type.")

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
	if (*binCells > 0 || *determineBins != "") && *cellTypeAnalysis != "" {
		log.Fatalf("Bin cells cannot be used with -cellTypeAnalysis")
	}
	if *altMapping != "" && *samOut != "" {
		log.Fatalf("altMapping is not compatable with samOut. If you are interested in collapsing UMI in a sam file use samFilter -collapseUMI")
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
		CountMatrix:          *countMatrix,
		NoOut:                *noOut,
		AltMapping:           *altMapping,
		CountMatrixCellTypes: *countMatrixCellTypes,
	}

	//output or input sequencing
	if *inputSequencing != "" {
		starrSeq.ParseInputSequencingSam(s)
	} else if *altMapping != "" {
		starrSeq.Alt(s)
	} else {
		parseCellrangerBam(s)
	}

}
