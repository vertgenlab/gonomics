package main

import (
	"log"
	"flag"
	"fmt"
	"strings"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
)

func vcfToBed(infile string, outfile string, delimiter string) {
	ch := vcf.GoReadToChan(infile)
	var note string
	var err error
	out := fileio.EasyCreate(outfile)
	defer out.Close()

	for v := range ch {
		words := strings.Split(v.Notes, "\t")
		note = words[0]
		for i := 0; i < len(words); i++ {
			note = note + delimiter + words[i]
		}
		_, err = fmt.Fprintf(out, "%s\t%v\t%v\t%s%s%s%s%s%s%v%s%s%s%s%s%s%s%s\n", v.Chr, v.Pos-1, v.Pos, v.Id, delimiter, v.Ref, delimiter, v.Alt, delimiter, v.Qual, delimiter, v.Filter, delimiter, v.Info, delimiter, v.Format, delimiter, note)
		common.ExitIfError(err)
	}
}

func usage() {
	fmt.Print(
	"vcfToBed: Converts vcf to bed format. Intended for subsequent use for vcf LiftOver.\n" +
	"Usage:\n" +
	"bedFilter input.vcf output.bed\n" +
	"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var delimiter *string = flag.String("delimiter", "&", "Sets the output name column delimiter.")
	
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

	vcfToBed(infile, outfile, *delimiter)
}
