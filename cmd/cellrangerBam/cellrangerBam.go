package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func parseBam(inSam string, outTable string) {
	var k int = 0
	var bit int32
	var constructName, cellString string
	var output []string

	ch, _ := sam.GoReadToChan(inSam)
	out := fileio.EasyCreate(outTable)

	for i := range ch {
		num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags
		bit = num.(int32)

		if bit&8 == 8 { // bit 8 is the flag for UMI used in final count.
			construct, _, _ := sam.QueryTag(i, "GX") // gene associated with UMI
			constructName = construct.(string)
			cell, _, _ := sam.QueryTag(i, "CB") // cell associated with UMI
			cellString = cell.(string)
			write := fmt.Sprintf("%s\t%s", cellString, constructName)
			fileio.WriteToFileHandle(out, write)
			k++
		}
	}
	fmt.Println("Found this many valid UMIs: ", k)
	
	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("cellrangerBam -- Takes in a cellranger bam file of STARR-seq reads and parses the extra flag field to pull out the" +
		"representitive read for each UMI and which construct it belongs to. It will also report the cell which it was found in.\n" +
		"Usage: \n" +
		"cellrangerBam inFile outFile\n\n")
	flag.PrintDefaults()
}
func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	a := flag.Arg(0)
	b := flag.Arg(1)
	parseBam(a, b)
}
