package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math/rand"
	"os"
	"sort"
	"strings"
)

func combineBams(a string, out *fileio.EasyWriter, iteration *int, length int, bw *sam.BamWriter) *sam.BamWriter {
	var err error
	var columns, columns2, columns3 []string
	var j string
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
				columns3 = strings.Split(j, "-")
				columns3[1] = fileio.IntToString(*iteration)
				join := strings.Join(columns3, "-")
				columns[n] = join
				join3 := strings.Join(columns, "\t")
				i.Extra = join3
			}
		}
		sam.WriteToBamFileHandle(bw, i, 0)
	}
	exception.PanicOnErr(err)
	*iteration++
	return bw
}

func binnedPseudobulk(inSlices [][]string, out *fileio.EasyWriter, norm string) {
	var j, i string
	var count float64
	var found bool
	var columns []string

	whichBin := 'A'

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

	binnedCells := make([][]string, numBins) //make a slicce of slice of string with the size of the user-specified number of bins
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
func singleCellAnalysis(cellTypeSlice []string, cellTypeInfo string, allConstructs []string, normalize string, writer *fileio.EasyWriter) {
	var found, found2 bool
	var cellType, j, m, cluster string
	var count float64
	var outMatrix, columns, columns2 []string

	cellTypeMap := make(map[string]string) // [cellBarcode]cellType
	allCellTypes := make(map[string]int)   // [cellType]placeholderInt

	cto := fileio.Read(cellTypeInfo) // read in the tab delimited file with cell and cell type info
	for _, i := range cto {          //make a map where the cell barcode is the key and the cell type is the value
		columns = strings.Split(i, "\t")
		cellTypeMap[columns[0]] = columns[1]
		_, found2 = allCellTypes[columns[1]] // also populate a map that contains a list of all cell types
		if !found2 {
			allCellTypes[columns[1]] = 1
		}
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
		if normalize != "" { // normalize the singleCellCount map if needed
			inputNormalize(singleCellTypeMap, normalize)
		}
		for m = range singleCellTypeMap { // write out the results from that cell type
			write := fmt.Sprintf("%s\t%s\t%f", cellType, m, singleCellTypeMap[m]) //cell type \t construct \t count
			outMatrix = append(outMatrix, write)
		}
		fmt.Println(fmt.Sprintf("Found %d raw counts in the cell type: %s", a, cellType))
	}
	sort.Strings(outMatrix)       //sort the slice before writing
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
	if len(inputNormValues) != len(mp) {
		fmt.Println("The input normalization table doesn't have the same number of constructs as was found in the input bam.")
		// trying to find the best way to throw an error if there is a raw count value that doesn't get normalized
	}
	for _, i := range inputNormValues { //iterate over the table with input normalization values and use those to edit the counts in the map
		columns = strings.Split(i, "\t")
		total, found = mp[columns[0]]
		if found {
			mp[columns[0]] = total * parse.StringToFloat64(columns[1])
		} else {
			fmt.Println("Construct not found in map for normalization: ", columns[0])
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
	for _, i := range writeSlice {
		fileio.WriteToFileHandle(writer, i)
	}
}

// parseBam takes in a cellranger bam file and pulls out reads that are representative of the UMI and also returns the construct associated with the UMI.
func parseBam(inSam string, outTable string, byCell bool, normalize string, samOut bool, cellTypeAnalysis string, numBins int) {
	var k int = 0
	var constructName, cellString, cellByConstructName string
	var count float64
	var found bool
	var bit int32
	var norm bool = false
	var sc bool = false
	var noSettings bool = false
	var allConstructs, cellTypeSlice []string

	if normalize != "" { //create a bool for normalize
		norm = true
	}
	if cellTypeAnalysis != "" { //create a bool for cellTypeAnalysis
		sc = true
	}
	if !sc && numBins < 1 && !samOut && !byCell {
		noSettings = true
	}

	ch, head := sam.GoReadToChan(inSam)

	out := fileio.EasyCreate(outTable)

	if samOut {
		sam.WriteHeaderToFileHandle(out, head)
	}

	pseudobulkMap := make(map[string]float64)

	for i := range ch { //interate of over the chanel of sam.Sam
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
		bit = num.(int32)
		if bit&8 == 8 { // bit 8 is the flag for a UMI that was used in final count (I call these "valid" UMIs.
			k++
			construct, _, _ := sam.QueryTag(i, "GX") // get the construct associated with valid UMI
			constructName = construct.(string)
			if byCell || sc || numBins > 0 { //byCell or singleCellAnalysis or binCells
				cell, _, _ := sam.QueryTag(i, "CB") // get the cell barcode associated with valid UMI
				cellString = cell.(string)
				cellByConstructName = fmt.Sprintf("%s\t%s", cellString, constructName) //list of all construct reads and what cell barcodes they belong to
				if byCell {                                                            //if you only want to print out the construct and the cell it was found in
					fileio.WriteToFileHandle(out, cellByConstructName)
				}
				if sc || numBins > 0 { //get some data for the singleCellAnalysis and distributeCells functions
					allConstructs = append(allConstructs, constructName) //list of all constructs found in the bam
					cellTypeSlice = append(cellTypeSlice, cellByConstructName)
				}
				continue
			} else if samOut { //write output as sam
				sam.WriteToFileHandle(out, i)
				continue
			} else if numBins > 0 {
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
	fmt.Println("Found this many valid UMIs: ", k)
	if sc { //singleCellAnalysis, handles all the writing as well
		singleCellAnalysis(cellTypeSlice, cellTypeAnalysis, allConstructs, normalize, out)
	}
	if numBins > 0 { //goes to binning and pseudobulk. Handles normalization and writing as well
		binnedCells := distributeCells(cellTypeSlice, numBins)
		binnedPseudobulk(binnedCells, out, normalize)
	}
	if norm && noSettings { //normalize pseudobulk
		inputNormalize(pseudobulkMap, normalize)
	}
	if noSettings { //write out pseudobulk
		writeMap(pseudobulkMap, out)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("cellrangerBam -- Takes in a cellranger bam file of STARR-seq reads and parses the extra flags field to pull out the" +
		"representative read for each UMI and which construct it belongs to. The output is a tab-delimited table of read-counts for each constructs.\n" +
		"NOTE: This function works best with STARR-seq libraries where constructs don't have much similarity with each other.\n" +
		"For libraries that need barcoding (like GWAS or cross-species comparisons) it is best practice to use samFilter and scCount" +
		"with a GTF corresponding to construct barcodes. \n" +
		"Usage: \n" +
		"cellrangerBam [options] inFile outFile\n\n")
	flag.PrintDefaults()
}

func main() {
	var byCell *bool = flag.Bool("byCell", false, "Will report the construct that each UMI belongs to and which cell in which it was found in a tab-delimited table.")
	var normalize *string = flag.String("normalize", "", "Takes in a tab delimited table with construct name and input normalization value. Must be used with -pseudobulk.")
	var samOut *bool = flag.Bool("samOut", false, "Output will be the reads that have valid UMIs in sam format")
	var cellTypeAnalysis *string = flag.String("cellTypeAnalysis", "", "Takes in a tab delimited file that has cell barcode and cell type identification. The ouptut of options will be a matrix that has counts for each construct in each cell type. The Seurat command WhichCells() can be used to generate the required list.")
	var binCells *int = flag.Int("binCells", 0, "Number of bins to randomly assign cells to. The output will be a psudobulk table for each bin")
	var combineGEMs *string = flag.String("combineGEMs", "", "Combine multiple GEM wells from the same experiment STARR-seq experiment to be analysed together. Will accept a comma separated list of bam files. The header from the argument 1 bam file input will be the one used. NOTE: This should only be used as an input for this program, not before any Seurat or similar analyses.")

	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *byCell && (*normalize != "" || *samOut) {
		log.Fatalf("Error: byCell cannot be used with normalize or samOut.")
	}

	if *normalize != "" && *samOut {
		log.Fatalf("Error: normalize and samOut cannot be used together.")
	}

	if (*byCell || *samOut) && *cellTypeAnalysis != "" {
		log.Fatalf("Error: cellTypeAnalysis is incompatable with byCell and samOut.")
	}

	if *binCells < 0 {
		log.Fatalf("Error: -binCells must be a positive intiger")
	}

	if *cellTypeAnalysis != "" && *combineGEMs != "" {
		log.Fatalf("Error: cellTypeAnalysis and combineGEMs cannot be used together currently.")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	a := flag.Arg(0)
	b := flag.Arg(1)

	var bw *sam.BamWriter
	var combineBamWriter *fileio.EasyWriter
	var err error
	var path []string
	if *combineGEMs != "" {
		var tmpFileSlice []string
		tmpFileSlice = append(tmpFileSlice, a)
		bams := strings.Split(*combineGEMs, ",")
		for _, i := range bams {
			path = strings.Split(i, "/")
			tmpFileSlice = append(tmpFileSlice, path[len(path)-1])
		}
		tmpFileSlice = append(tmpFileSlice, ".bam")
		tmpFileName := strings.Join(tmpFileSlice, "_")
		combineBamWriter = fileio.EasyCreate(tmpFileName)
		p := 1
		bw = combineBams(a, combineBamWriter, &p, len(bams), bw)
		for _, i := range bams {
			combineBams(i, combineBamWriter, &p, len(bams), bw)
		}
		a = tmpFileName

		err = bw.Close()
		exception.PanicOnErr(err)
		err = combineBamWriter.Close()
		exception.PanicOnErr(err)

	}

	parseBam(a, b, *byCell, *normalize, *samOut, *cellTypeAnalysis, *binCells)

	if *combineGEMs != "" {
		err = os.Remove(a)
		exception.PanicOnErr(err)
	}
}
