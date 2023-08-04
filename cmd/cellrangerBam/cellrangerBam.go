package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sort"
	"strings"
)

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
	for _, i := range inputNormValues {
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

// cellrangerBam takes in a cellranger bam file and pulls out reads that are representative of the UMI and also returns the construct associated with the UMI.
func cellrangerBam(inSam string, outTable string, byCell bool, normalize string, samOut bool) {
	var k int = 0
	var constructName, cellString, write string
	var count float64
	var found bool
	var bit uint8
	var norm bool = false

	if normalize != "" {
		norm = true
	}

	ch, head := sam.GoReadToChan(inSam)

	out := fileio.EasyCreate(outTable)

	if samOut {
		sam.WriteHeaderToFileHandle(out, head)
	}
	mp := make(map[string]float64)

	for i := range ch {
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags
		bit = num.(uint8)
		if bit&8 == 8 { // bit 8 is the flag for a UMI that was used in final count.
			k++
			construct, _, _ := sam.QueryTag(i, "GX") // gene associated with UMI
			constructName = construct.(string)
			if byCell { //byCell
				cell, _, _ := sam.QueryTag(i, "CB") // cell associated with UMI
				cellString = cell.(string)
				write = fmt.Sprintf("%s\t%s", cellString, constructName)
				fileio.WriteToFileHandle(out, write)
				continue
			} else if samOut { //write output as sam
				sam.WriteToFileHandle(out, i)
				continue
			}
			//pseudobulk, default behavior
			count, found = mp[constructName]
			if !found {
				mp[constructName] = 1
			} else {
				mp[constructName] = count + 1
			}
		}
	}
	fmt.Println("Found this many valid UMIs: ", k)
	if norm {
		inputNormalize(mp, normalize)
	}
	writeMap(mp, out)

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

	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *byCell && *normalize != "" {
		log.Fatalf("Error: byCell and normalize cannot be used together.")
	}

	if *byCell && *samOut {
		log.Fatalf("Error: byCell and samOut cannot be used together.")
	}

	if *normalize != "" && *samOut {
		log.Fatalf("Error: normalize and samOut cannot be used together.")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	a := flag.Arg(0)
	b := flag.Arg(1)

	cellrangerBam(a, b, *byCell, *normalize, *samOut)
}
