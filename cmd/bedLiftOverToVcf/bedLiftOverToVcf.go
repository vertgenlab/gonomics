// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

func bedLiftOverToVcf(infile string, outfile string, delimiter string) {
	ch := bed.GoReadToChan(infile)
	var err error
	out := fileio.EasyCreate(outfile)
	defer out.Close()

	for v := range ch {
		var output string = fmt.Sprintf("%s\t%v", v.Chrom, v.ChromEnd)
		words := strings.Split(v.Name, delimiter)
		for i := 0; i < len(words); i++ {
			output = output + "\t" + words[i]
		}
		_, err = fmt.Fprintf(out, output+"\n")
		common.ExitIfError(err)
	}

}

func usage() {
	fmt.Print(
		"bedLiftOverToVcf: Converts a bed liftOver product back to a vcf. Intended for vcf Liftover.\n" +
			"Usage:\n" +
			"bedFilter input.bed output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var delimiter *string = flag.String("delimiter", "&", "Sets the input name column delimiter.")

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

	bedLiftOverToVcf(infile, outfile, *delimiter)
}
