package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math/rand"
	"os"
	"sort"
	"strings"
)

type Settings struct {
	inFile         string
	outFile        string
	inputNormalize string
	byCell         bool
	samOut         bool
	scAnalysis     string
	binCells       int
	umiSat         bool
	gfpNorm        string
	bed            string
	ncNorm         string
}

//getNcAvg takes in the input normalized matrix of single cell data, the list of negative control sequences and an empty slice of slices and will return the average negative control reads for each
//cell type
func getNcAvg(matrix []string, ncSlice []string, ncVals [][]float64) []float64 {
	var found bool
	var currCluster string
	var columns []string
	var prevCluster string
	var ncAvg []float64

	columns = strings.Split(matrix[0], "\t")
	currCluster = columns[0]

	clusterIdx := 0
	for p, i := range matrix {
		found = false
		columns = strings.Split(i, "\t")
		currCluster = columns[0]
		if p > 0 && currCluster != prevCluster {
			clusterIdx++
		}
		for _, j := range ncSlice {
			if j == columns[1] {
				found = true
				break
			}
		}
		if found {
			ncVals[clusterIdx] = append(ncVals[clusterIdx], parse.StringToFloat64(columns[2]))
		}
		prevCluster = currCluster
	}
	for i, _ := range ncVals {
		ncAvg = append(ncAvg, numbers.AverageFloat64(ncVals[i]))
	}
	return ncAvg
}

//normScToNegativeCtrls is the initial function to normalize single-cell data to negative control reads in the same cell type. It takes in a matrix that is the input-normalized, the settings struct,
//the number of cell types and returns the negative control normalized read counts
func normScToNegativeCtrls(matrix []string, s Settings, numCellTypes int) []string {
	var ncSlice, columns []string
	var prevCluster, newLine string

	ncVals := make([][]float64, numCellTypes)
	nc := fileio.Read(s.ncNorm)
	for _, i := range nc {
		ncSlice = append(ncSlice, i)
	}
	ncAvg := getNcAvg(matrix, ncSlice, ncVals)

	columns = strings.Split(matrix[0], "\t")
	currCluster := columns[0]

	clusterIdx := 0
	for j, i := range matrix {
		columns = strings.Split(i, "\t")
		currCluster = columns[0]
		if j > 0 && currCluster != prevCluster {
			clusterIdx++
		}
		columns[2] = fmt.Sprintf("%f", parse.StringToFloat64(columns[2])/ncAvg[clusterIdx])
		newLine = strings.Join(columns, "\t")
		matrix[j] = newLine
		prevCluster = currCluster
	}
	return matrix
}

//gfpNormalize iterates over a single-cell counts map and applies a GFP normalization factor.
func gfpNormalize(countsMap map[string]float64, normFactor float64) {
	var counts float64
	for i := range countsMap {
		counts, _ = countsMap[i]
		countsMap[i] = counts * normFactor
	}
}

//gfpNormFactor takes in a map of cellType-gfpReads and returns a map with cellType-gfpNormalizationValue. Normalization value is calculated by the percentage of GFP reads if all
//cell types were equal, and dividing it by the observed percentage of GFP reads in a given cell type
func gfpNormFactor(clusterGFP map[string]int) map[string]float64 {
	var totalCounts, clusterCounts int
	var actualPerc float64
	gfpNormMap := make(map[string]float64)
	idealPerc := 1.0 / float64(len(clusterGFP))
	for i := range clusterGFP {
		clusterCounts, _ = clusterGFP[i]
		totalCounts = clusterCounts + totalCounts
	}
	for i := range clusterGFP {
		actualPerc = float64(clusterGFP[i]) / float64(totalCounts)
		gfpNormMap[i] = idealPerc / actualPerc
	}
	return gfpNormMap
}

//parseGfpBam is similar to the parseBam function but handles bams containing GFP reads. This function also takes in a map of cellBx-cellType and an empty map containing each cellType
//It returns a map that contains [cellType] = gfpReads
func parseGfpBam(gfpBam string, cellTypeMap map[string]string, clusterGFP map[string]int) map[string]int {
	var bit int32
	var cluster string
	var count int
	var gfpStats []string
	var found bool

	inChan, _ := sam.GoReadToChan(gfpBam)
	for i := range inChan {
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
		bit = num.(int32)
		if bit&8 == 8 {
			cellBx, _, _ := sam.QueryTag(i, "CB")
			cluster, found = cellTypeMap[cellBx.(string)]
			if found {
				count, _ = clusterGFP[cluster]
				clusterGFP[cluster] = count + 1
			}
		}
	}
	//print out some GFP read counts / cluster
	var gfpReads int
	for i := range clusterGFP {
		gfpReads, _ = clusterGFP[i]
		gfpStats = append(gfpStats, fmt.Sprintf("found %d GFP reads in %s", gfpReads, i))
	}
	for _, i := range gfpStats {
		fmt.Println(i)
	}
	return clusterGFP
}

// umiSaturation randomly subsets the whole bam file (10% to 100% of all reads) and calculates how many UMIs are in those subests. The output is a tab delimited text file.
func umiSaturation(umiBxSlice []string, out *fileio.EasyWriter) {
	var perc, randNum float64
	var j string
	var count int

	fileio.WriteToFileHandle(out, "totalReads\tumis")

	for i := 1; i <= 10; i++ {
		perc = float64(i) / 10
		mp := make(map[string]bool)
		count = 0
		for _, j = range umiBxSlice {
			randNum = rand.Float64()
			if randNum > perc {
				continue
			}
			count++
			mp[j] = true
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%d", count, len(mp)))
	}
}

// combineBams is triggered when there are more than 1 input bams in a comma delimited list. This function appends a "_int" to the back of the cell barcode correstponding to the input file order
// and then concatenates all bam files together
func combineBams(a string, out *fileio.EasyWriter, iteration *int, length int, bw *sam.BamWriter) *sam.BamWriter {
	var err error
	var columns, columns2 []string
	var j, join string
	var n int

	inChan, head := sam.GoReadToChan(a)
	if *iteration == 1 {
		bw = sam.NewBamWriter(out, head)
	}
	for i := range inChan {
		err = sam.ParseExtra(&i)
		columns = strings.Split(i.Extra, "\t")
		for n, j = range columns {
			columns2 = strings.Split(j, ":")
			if columns2[0] == "CB" {
				columns[n] = fmt.Sprintf("%s_%d", j, *iteration)
				join = strings.Join(columns, "\t")
				i.Extra = join
			}
		}
		sam.WriteToBamFileHandle(bw, i, 0)
	}
	exception.PanicOnErr(err)
	*iteration++
	return bw
}

// binnedPseudobulk is the psuedobulk function if the -binCells option is used. It adds an addition column to the dataframe corresponding to bin identity
func binnedPseudobulk(inSlices [][]string, out *fileio.EasyWriter, norm string) {
	var j, i string
	var count float64
	var found bool
	var columns []string

	whichBin := 'A'
	fileio.WriteToFileHandle(out, "construct\tcounts\tbin")
	for _, bin := range inSlices {
		mp := make(map[string]float64)
		var toWrite []string
		for _, j = range bin {
			columns = strings.Split(j, "\t")
			count, found = mp[columns[1]]
			if !found {
				mp[columns[1]] = 1
			} else {
				mp[columns[1]] = count + 1
			}
		}
		if norm != "" {
			inputNormalize(mp, norm)
		}
		for i = range mp {
			toWrite = append(toWrite, fmt.Sprintf("%s\t%f\t%c", i, mp[i], whichBin))
		}
		sort.Strings(toWrite)
		for _, i = range toWrite {
			fileio.WriteToFileHandle(out, i)
		}
		whichBin++
	}
}

// determineBin takes in an int and the probablity of any given bin and determines the bin the cell belongs to.
func determineBin(numBins int, prob float64) int {
	var bin int

	randNum := rand.Float64()      // draw a random number
	for j := 0; j < numBins; j++ { //iterate over the number of bins
		if j == numBins { //if we've gotten to the last bin and the cell still hasn't been partitioned, we don't have to check anything else we can just put it there.
			return j
		}
		if randNum >= prob*float64(j) && randNum < prob*float64(j+1) { // check to see if our random number belongs to that bin
			bin = j
		}
	}
	return bin
}

// distributeCells takes in a slice of string that is created in parseBam that has a list of cellBarcodes and constructs. It will partition those cells and constructs into separate slices. The number of bins is a user-input variable
func distributeCells(cellTypeSlice []string, numBins int) [][]string {
	var currCell string
	var columns []string
	var bin, count int
	var found bool
	whichBin := 'A'                   //for counting cells per bin
	whichBinMap := make(map[rune]int) //for counting cells per bin

	binnedCells := make([][]string, numBins) //make a slice of slice of string with the size of the user-specified number of bins
	prob := 1.0 / float64(numBins)           // determine the probablity that the cell belongs to a particular bin

	sort.Strings(cellTypeSlice) //sort the slice of strings containing (cellBarcode \t construct) so that indentical cell barcodes line up next to one annother

	for _, i := range cellTypeSlice {
		columns = strings.Split(i, "\t")
		if columns[0] == currCell { // if the cell barcode is the same as the one the loop just saw, partition that cell-construct pair into the same bin as the one before
			binnedCells[bin] = append(binnedCells[bin], i)
			continue
		}
		bin = determineBin(numBins, prob)              //function to determine which bin to put the new cell into.
		binnedCells[bin] = append(binnedCells[bin], i) //put cell into correct bin
		currCell = columns[0]                          //set the cell barcode that was just partitioned to current cell for the next iteration
		// next part of the code is just to count how many cells went into each bin and print out those values. May keep it in, but it is also a check to make sure things are running correctly
		count, found = whichBinMap[whichBin+rune(bin)]
		if !found {
			whichBinMap[whichBin+rune(bin)] = 1
		} else {
			whichBinMap[whichBin+rune(bin)] = count + 1
		}
	}
	var outSlice []string
	for i := range whichBinMap {
		outSlice = append(outSlice, fmt.Sprintf("Bin: %c\tcells: %d", i, whichBinMap[i]))
	}
	sort.Strings(outSlice)
	for _, i := range outSlice {
		fmt.Println(i)
	}
	return binnedCells
}

// singleCellAnalysis will create a count for each construct in each cell type. This function takes
func singleCellAnalysis(s Settings, cellTypeSlice []string, allConstructs []string, writer *fileio.EasyWriter) {
	var found, found2 bool
	var cellType, j, m, cluster string
	var count float64
	var outMatrix, columns, columns2 []string
	var gfpNormFactorMap map[string]float64

	cellTypeMap := make(map[string]string) // [cellBarcode]cellType
	allCellTypes := make(map[string]int)   // [cellType]placeholderInt
	//countsPerCellType := make(map[string]int)

	cto := fileio.Read(s.scAnalysis) // read in the tab delimited file with cell and cell type info
	for _, i := range cto {          //make a map where the cell barcode is the key and the cell type is the value
		columns = strings.Split(i, "\t")
		cellTypeMap[columns[0]] = columns[1]
		_, found2 = allCellTypes[columns[1]] // also populate a map that contains a list of all cell types
		if !found2 {
			allCellTypes[columns[1]] = 0
		}
	}

	if s.gfpNorm != "" {
		gfpClusterMap := parseGfpBam(s.gfpNorm, cellTypeMap, allCellTypes)
		gfpNormFactorMap = gfpNormFactor(gfpClusterMap)
	}

	for cellType = range allCellTypes { //loop over the list of cell types
		var a int = 0
		singleCellTypeMap := make(map[string]float64) // make a map which will hold the construct counts for the current cell type ([construct]counts)
		for _, j = range allConstructs {              // loop over all construct that were found in the bam file to make every construct start at 0 reads
			singleCellTypeMap[j] = 0
		}
		for _, i := range cellTypeSlice { //loop over the slice created in cellrangerBam() that that contains: cellBarcode \t construct
			columns2 = strings.Split(i, "\t")
			cluster, found = cellTypeMap[columns2[0]] //search for that cell barcode in the cellType map
			if found && cluster == cellType {         // if that cell is found and the cell type in the map matches the cell type we are looping through, search for the construct corresponding with that cellBarcode in the singleCellCount map and add one to the value
				count, _ = singleCellTypeMap[columns2[1]]
				singleCellTypeMap[columns2[1]] = count + 1
				a++
			}
		}
		if s.inputNormalize != "" { // normalize the singleCellCount map if needed
			inputNormalize(singleCellTypeMap, s.inputNormalize)
		}
		if s.gfpNorm != "" {
			gfpNormalizationFactor := gfpNormFactorMap[cellType]
			gfpNormalize(singleCellTypeMap, gfpNormalizationFactor)
		}
		for m = range singleCellTypeMap { // write out the results from that cell type
			write := fmt.Sprintf("%s\t%s\t%f", cellType, m, singleCellTypeMap[m]) //cell type \t construct \t count
			outMatrix = append(outMatrix, write)
		}
		fmt.Println(fmt.Sprintf("Found %d raw counts in the cell type: %s", a, cellType))
	}
	sort.Strings(outMatrix) //sort the slice before writing

	if s.ncNorm != "" {
		outMatrix = normScToNegativeCtrls(outMatrix, s, len(allCellTypes))
	}
	if s.ncNorm != "" {
		fileio.WriteToFileHandle(writer, "cellCluster\tconstruct\tcounts/ncCounts")
	} else {
		fileio.WriteToFileHandle(writer, "cellCluster\tconstruct\tcounts")
	}
	for _, i := range outMatrix { //once all cell types have been looped through, write out the total count matrix
		fileio.WriteToFileHandle(writer, i)
	}
}

// inputNormalize takes in the psuedobulk map and an input normalization table and normalizes all the raw count values in the map
func inputNormalize(mp map[string]float64, normalize string) {
	var total float64
	var found bool
	var columns []string

	inputNormValues := fileio.Read(normalize)
	if len(inputNormValues) < len(mp) {
		fmt.Println("The input normalization table has less constructs than were found in the input bam. 1 or more constructs won't be normalized. Please check your input files.")
	} else if len(inputNormValues) > len(mp) {
		fmt.Println("The input normalization table has more constructs than were found in the input bam. Constructs not found in the bam will have 0 counts in the output.")
	}
	for _, i := range inputNormValues { //iterate over the table with input normalization values and use those to edit the counts in the map
		columns = strings.Split(i, "\t")
		total, found = mp[columns[0]]
		if found {
			mp[columns[0]] = total * parse.StringToFloat64(columns[1])
		} else {
			mp[columns[0]] = 0.0 //user wants to normalize this construct meaning that it was in the library. however, no reads were recovered for that construct so we add it to the map and set the value to zero
		}
	}
}

// writeMap simply writes out the pseudobulk and/or input normalized values to an io.writer
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
func parseBam(s Settings) {
	var k int = 0
	var constructName, cellString, cellByConstructName, umiBx, constructString string
	var count float64
	var found bool
	var bit int32
	var norm bool = false
	var sc bool = false
	var noSettings bool = false
	var allConstructs, cellTypeSlice, umiBxSlice, constructSlice []string

	if s.inputNormalize != "" { //create a bool for normalize
		norm = true
	}
	if s.scAnalysis != "" { //create a bool for cellTypeAnalysis
		sc = true
	}
	if !sc && s.binCells < 1 && !s.samOut && !s.byCell && !s.umiSat {
		noSettings = true
	}

	ch, head := sam.GoReadToChan(s.inFile)

	out := fileio.EasyCreate(s.outFile)

	if s.samOut {
		sam.WriteHeaderToFileHandle(out, head)
	}

	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	if s.bed != "" {
		bedEntries := bed.Read(s.bed)
		for _, i := range bedEntries {
			tree = append(tree, i)
		}
		selectTree = interval.BuildTree(tree)
	}

	pseudobulkMap := make(map[string]float64)

	for i := range ch { //iterate of over the chanel of sam.Sam
		if s.umiSat {
			umi, _, _ := sam.QueryTag(i, "UB")
			cb, _, _ := sam.QueryTag(i, "CB")
			umiBx = fmt.Sprintf("%s_%s", umi, cb)
			umiBxSlice = append(umiBxSlice, umiBx)
		}
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
		bit = num.(int32)
		if bit&8 == 8 { // bit 8 is the flag for a UMI that was used in final count. I call these "valid" UMIs.
			k++
			if s.bed != "" {
				overlappedBedEntry := interval.Query(selectTree, i, "any")
				if len(overlappedBedEntry) != 1 {
					continue
				}
				constructString = fmt.Sprint(overlappedBedEntry[0])
				constructSlice = strings.Split(constructString, "\t")
				constructName = constructSlice[3]
			} else {
				construct, _, _ := sam.QueryTag(i, "GX") // get the construct associated with valid UMI
				constructName = construct.(string)
			}
			if s.byCell || sc || s.binCells > 0 { //byCell or singleCellAnalysis or binCells
				cell, _, _ := sam.QueryTag(i, "CB") // get the cell barcode associated with valid UMI
				cellString = cell.(string)
				cellByConstructName = fmt.Sprintf("%s\t%s", cellString, constructName) //list of all construct reads and what cell barcodes they belong to
				if s.byCell {                                                          //if you only want to print out the construct and the cell it was found in
					fileio.WriteToFileHandle(out, cellByConstructName)
				}
				if sc || s.binCells > 0 { //get some data for the singleCellAnalysis and distributeCells functions
					allConstructs = append(allConstructs, constructName) //list of all constructs found in the bam
					cellTypeSlice = append(cellTypeSlice, cellByConstructName)
				}
				continue
			} else if s.samOut { //write output as sam
				sam.WriteToFileHandle(out, i)
				continue
			} else if s.binCells > 0 {
				continue
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
	if s.umiSat {
		umiSaturation(umiBxSlice, out)
	}
	fmt.Println("Found this many valid UMIs: ", k)
	if sc { //singleCellAnalysis, handles all the writing as well
		singleCellAnalysis(s, cellTypeSlice, allConstructs, out)
	}
	if s.binCells > 0 { //goes to binning and pseudobulk. Handles normalization and writing as well
		binnedCells := distributeCells(cellTypeSlice, s.binCells)
		binnedPseudobulk(binnedCells, out, s.inputNormalize)
	}
	if norm && noSettings { //normalize pseudobulk
		inputNormalize(pseudobulkMap, s.inputNormalize)
	}
	if noSettings { //write out pseudobulk
		writeMap(pseudobulkMap, out)
	}
	err := out.Close()
	exception.PanicOnErr(err)
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
	var byCell *bool = flag.Bool("byCell", false, "Will report the construct that each UMI belongs to and which cell in which it was found in a tab-delimited table.")
	var inputNorm *string = flag.String("inputNorm", "", "Takes in a tab delimited table with construct name and input normalization value")
	var samOut *bool = flag.Bool("samOut", false, "Output will be the reads that have valid UMIs in sam format")
	var cellTypeAnalysis *string = flag.String("cellTypeAnalysis", "", "Takes in a tab delimited file that has cell barcode and cell type identification. "+
		"The ouptut of options will be a matrix that has counts for each construct in each cell type. The Seurat command WhichCells() can be used to generate the required list.")
	var binCells *int = flag.Int("binCells", 0, "Number of bins to randomly assign cells to. The output will be a psudobulk table for each bin")
	var umiSat *bool = flag.Bool("umiSat", false, "Create a UMI saturation curve of all reads in the input bam file.")
	var gfpNorm *string = flag.String("gfpNorm", "", "Bam file from the same scSTARR-seq experiment containing cellranger count alignments to GFP for cell cluster GFP normalization. "+
		"Multiple bam files from different GEM wells can be provided in a comma-separated list")
	var bed *string = flag.String("bed", "", "Use a bed file for assigning reads to constructs instead of the GTF that was used in cellranger mkref. The bed file must have the "+
		"name of the construct in the fourth field. Recammended for constructs with barcodes.")
	var ncNorm *string = flag.String("ncNorm", "", "Reports single cell data as the ratio the reads for each construct to the negative control reads in the same cell type."+
		"A file that conatins a line-delimited list of negative controls must be provided.")

	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *ncNorm != "" && *gfpNorm != "" {
		log.Fatalf("Error: only one normalization method for single-cell data must be used")
	}

	if (*ncNorm != "" || *gfpNorm != "") && *cellTypeAnalysis == "" {
		log.Fatalf("Error: -ncNorm or -gfpNorm must be used with -cellTypeAnalysis")
	}

	if *byCell && (*inputNorm != "" || *samOut) {
		log.Fatalf("Error: byCell cannot be used with normalize or samOut.")
	}

	if *inputNorm != "" && *samOut {
		log.Fatalf("Error: normalize and samOut cannot be used together.")
	}

	if (*byCell || *samOut) && *cellTypeAnalysis != "" {
		log.Fatalf("Error: cellTypeAnalysis is incompatable with byCell and samOut.")
	}

	if *binCells < 0 {
		log.Fatalf("Error: -binCells must be a positive intiger")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	var s Settings = Settings{
		inFile:         flag.Arg(0),
		outFile:        flag.Arg(1),
		inputNormalize: *inputNorm,
		byCell:         *byCell,
		scAnalysis:     *cellTypeAnalysis,
		binCells:       *binCells,
		umiSat:         *umiSat,
		samOut:         *samOut,
		gfpNorm:        *gfpNorm,
		bed:            *bed,
		ncNorm:         *ncNorm,
	}

	var bw *sam.BamWriter
	var combineBamWriter *fileio.EasyWriter
	var err error
	inFiles := strings.Split(s.inFile, ",")
	if len(inFiles) > 1 {
		var path, tmpFileSlice []string
		if *cellTypeAnalysis != "" {
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -cellTypeAnalysis. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		}
		for _, i := range inFiles {
			path = strings.Split(i, "/")
			tmpFileSlice = append(tmpFileSlice, path[len(path)-1])
		}
		tmpFileSlice = append(tmpFileSlice, ".bam")
		tmpFileName := strings.Join(tmpFileSlice, "_")
		combineBamWriter = fileio.EasyCreate(tmpFileName)
		p := 1
		for _, i := range inFiles {
			bw = combineBams(i, combineBamWriter, &p, len(inFiles), bw)
		}
		s.inFile = tmpFileName
		err = bw.Close()
		exception.PanicOnErr(err)
		err = combineBamWriter.Close()
		exception.PanicOnErr(err)

	}

	var gfpBams []string
	if *gfpNorm != "" {
		var bwGFP *sam.BamWriter
		var combineBamWriterGFP *fileio.EasyWriter
		gfpBams = strings.Split(*gfpNorm, ",")
		if len(gfpBams) > 1 {
			var path, tmpFileSlice []string
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -gfpNorm. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
			for _, i := range gfpBams {
				path = strings.Split(i, "/")
				tmpFileSlice = append(tmpFileSlice, path[len(path)-1])
			}
			tmpFileSlice = append(tmpFileSlice, ".bam")
			tmpFileName := strings.Join(tmpFileSlice, "_")
			combineBamWriterGFP = fileio.EasyCreate(tmpFileName)
			p := 1
			for _, i := range gfpBams {
				bwGFP = combineBams(i, combineBamWriterGFP, &p, len(gfpBams), bwGFP)
			}
			s.gfpNorm = tmpFileName
			err = bwGFP.Close()
			exception.PanicOnErr(err)
			err = combineBamWriterGFP.Close()
			exception.PanicOnErr(err)

		}
	}

	parseBam(s)

	if len(inFiles) > 1 {
		err = os.Remove(s.inFile)
		exception.PanicOnErr(err)
	}
	if len(gfpBams) > 1 {
		err = os.Remove(s.gfpNorm)
		exception.PanicOnErr(err)
	}
}
