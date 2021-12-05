// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func afsPmf(vcfFile string, outFile string, UnPolarized bool) {
	g, err := popgen.VcfToAfs(vcfFile, UnPolarized)
	exception.FatalOnErr(err)
	out := fileio.EasyCreate(outFile)

	numAlleles := g.Sites[0].N
	frequencySlice := make([]int, numAlleles)

	for _, site := range g.Sites {
		if site.N == numAlleles {
			frequencySlice[site.I]++
		}
	}
	_, err = fmt.Fprintf(out, "AlleleFrequency\tNumAlleles\tFloatFrequency\n")
	exception.PanicOnErr(err)
	for i := 1; i < len(frequencySlice); i++ {
		_, err = fmt.Fprintf(out, "%v\t%v\t%f\n", i, frequencySlice[i], float64(frequencySlice[i])/float64(len(g.Sites)))
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"afsPmf - Returns the probability mass function of the allele frequency spectrum in an input vcf file.\n All segregating sites must have the same value for N. Useful for graphing.\n" +
			"afsPmf in.vcf out.txt\n")
}

func main() {
	var unPolarized *bool = flag.Bool("unPolarized", false, "afsPmf creates polarized derived frequency spectra by default. When true, the cmd returns unpolarized site frequency spectra.")

	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	afsPmf(vcfFile, outFile, *unPolarized)
}
