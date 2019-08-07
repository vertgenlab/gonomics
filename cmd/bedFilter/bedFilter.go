package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
	"strings"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
)

func bedFilter(infile string, outfile string, threshold int64) {
	var outlist []*bed.Bed
	var line string
	var startNum, endNum int64
	var doneReading bool = false
	var current *bed.Bed
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		startNum = common.StringToInt64(words[1])
		endNum = common.StringToInt64(words[2])
		if common.StringToInt64(words[4]) >= threshold {
			current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, Name: words[1], Score: common.StringToInt64(words[4])}
			outlist = append(outlist, current)
		}
	
	}
	bed.Write(outfile, outlist, 5)
	/* Old version, memory expensive

	var records []*bed.Bed = bed.Read(infile)
	var outlist []*bed.Bed

	for i := 0; i < len(records); i++ {
		if records[i].Score >= *threshold {
			outlist = append(outlist, records[i])
		}
	}

	bed.Write(outfile, outlist, 5)
	*/
}

func usage() {
	fmt.Print(
		"bedFilter - removes bed entries below a specified threshold score.\n" +
		"Usage:\n" +
		"bedFilter input.bed output.bed\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var threshold *int64 = flag.Int64("threshold", 0, "Specifies the threshold value")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedFilter(infile, outfile, *threshold)
}
