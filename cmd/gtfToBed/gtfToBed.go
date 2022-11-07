// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"log"
	"strings"
)

func gtfToBed(fileName string, outFile string, tss bool, chromSizeFile string) {
	var err error
	var line string
	var nameString string
	var currBed bed.Bed
	var doneReading bool = false

	if tss {
		g := gtf.Read(fileName)
		sizes := chromInfo.ReadToMap(chromSizeFile)
		beds := gtf.GenesToTssBed(g, sizes)
		bed.Write(outFile, beds)
	} else {
		file := fileio.EasyOpen(fileName)
		out := fileio.EasyCreate(outFile)
		for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
			words := strings.Split(line, "\t")
			nameString = words[1] + ":" + words[2]
			for i := 5; i < len(words); i++ {
				nameString = nameString + ":" + words[i]
			}
			currBed = bed.Bed{Chrom: words[0], ChromStart: common.StringToInt(words[3]) - 1, ChromEnd: common.StringToInt(words[4]), Name: nameString, Score: 0, Strand: bed.Positive, FieldsInitialized: 6}
			if words[6] == "-" {
				currBed.Strand = bed.Negative
			}
			bed.WriteBed(out, currBed)
		}
		err = file.Close()
		exception.PanicOnErr(err)
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print(
		"gtfToBed - Converts a gtf file into bed format.\n" +
			"Usage:\n" +
			"gtfToBed in.gtf out.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var tss *bool = flag.Bool("tss", false, "Return a bed of tss positions annotated only with the geneName. Must provide chrom sizes file.")
	var chromSizeFile *string = flag.String("chromSizeFile", "", "Specifies the name of a chrom.sizes file.")
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	fileName := flag.Arg(0)
	outFile := flag.Arg(1)

	gtfToBed(fileName, outFile, *tss, *chromSizeFile)
}
