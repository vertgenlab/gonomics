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

type Settings struct {
	UnPolarized             bool
	PlotSelectionLikelihood string
	LeftBound               float64
	RightBound              float64
	NumberOfPoints          int
	IntegralError           float64
	DivergenceAscertainment bool
	D                       int //size of the ascertainment subset
}

func vcfAfs(vcfFile string, outFile string, s Settings) {
	var err error
	g, err := popgen.VcfToAfs(vcfFile, s.UnPolarized, false) //hardcoded DivergenceAscertainment to false
	exception.FatalOnErr(err)
	f := popgen.AfsToFrequency(*g)
	out := fileio.EasyCreate(outFile)
	for i := 0; i < len(f); i++ {
		_, err = fmt.Fprintf(out, "%f\n", f[i])
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)

	if s.PlotSelectionLikelihood != "" {
		popgen.PlotAfsLikelihood(*g, s.PlotSelectionLikelihood, s.LeftBound, s.RightBound, s.NumberOfPoints, s.IntegralError, s.DivergenceAscertainment, s.D)
	}
}

func usage() {
	fmt.Print(
		"vcfAfs - Returns allele frequency spectrum information in a text file for graphing.\n" +
			"vcfAfs - vcfFile.vcf outFile.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var unPolarized *bool = flag.Bool("unPolarized", false, "vcfAfs creates polarized derived frequency spectra by default. When true, the cmd returns unpolarized site frequency spectra.")
	var plotSelectionLikelihood *string = flag.String("plotSelectionLikelihood", "", "Plots the likelihood function for the selection parameter from this VCF to a specified output file.")
	var leftBound *float64 = flag.Float64("leftBound", -10.0, "Left bound for plotting of the likelihood function.")
	var rightBound *float64 = flag.Float64("rightBound", 10.0, "Right bound for plotting of the likelihood function.")
	var numberOfPoints *int = flag.Int("numberOfPoints", 100, "Set the number of points at which to evaluate the likelihood function.")
	var integralError *float64 = flag.Float64("integralError", 1e-5, "Sets the error threshold for integral calculations in the likelihood calculation.")
	var divergenceAscertainment *bool = flag.Bool("divergenceAscertainment", false, "Corrects for divergence-based ascertainment bias in the likelihood calculation.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		UnPolarized:             *unPolarized,
		PlotSelectionLikelihood: *plotSelectionLikelihood,
		LeftBound:               *leftBound,
		RightBound:              *rightBound,
		NumberOfPoints:          *numberOfPoints,
		IntegralError:           *integralError,
		DivergenceAscertainment: *divergenceAscertainment,
		D:                       1, //hardcoded for now
	}

	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	vcfAfs(vcfFile, outFile, s)
}
