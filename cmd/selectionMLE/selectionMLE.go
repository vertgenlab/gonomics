package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func selectionMLE(inFile string, outFile string, s popgen.MleSettings) {
	data, err := popgen.VcfToAfs(inFile, s.UnPolarized, s.DivergenceAscertainment) //VcfToAFS is written with polarized as the argument for clarity, so the bool is flipped here.
	exception.FatalOnErr(err)
	answer := popgen.SelectionMaximumLikelihoodEstimate(*data, s)
	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out,"#FILENAME\tMaximumLikelihood\n")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "%s\t%e\n", inFile, answer)
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"selectionMLE - Performs maximum likelihood estimation of selection on variants from an input VCF format file.\n" +
			"All bases are assumed to have the  same selection parameter, and the maximum of the likelihood function is found with\n" +
			"golden-segment search.\n" +
			"Usage:\n" +
			"selectionMLE input.vcf out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var leftBound *float64 = flag.Float64("leftBound", -10, "Set the leftBound for the MLE search space.")
	var rightBound *float64 = flag.Float64("rightBound", 10, "Set the right bound for the MLE search space.")
	var errorThreshold *float64 = flag.Float64("errorThreshold", 1e-5, "Set the desired level of accuracy in the maximum likelihood estimate.")
	var expectedNumArgs int = 2
	var unPolarized *bool = flag.Bool("unPolarized", false, "Disable the requirement for ancestor annotation and use unpolarized site frequency spectrum. Use with caution.")
	var divergenceAscertainment *bool = flag.Bool("divergenceAscertainment", false, "Make a divergence-based ascertainment correction.")
	var integralError *float64 = flag.Float64("integralError", 1e-7, "Set the error threshold for numerical integration.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 or 2 to reveal different levels of debug print statements to standard output. Currently none available.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	log.SetFlags(0)
	flag.Parse()

	options := popgen.MleSettings{
		Left: *leftBound,
		Right: *rightBound,
		Error: *errorThreshold,
		UnPolarized:             *unPolarized,
		DivergenceAscertainment: *divergenceAscertainment,
		D:                       1, //D is hardcoded as 1 for now. This represents the size of the ascertainment subset.
		IntegralError:           *integralError,
		Verbose:                 *verbose,
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	selectionMLE(vcfFile, outFile, options)
}

