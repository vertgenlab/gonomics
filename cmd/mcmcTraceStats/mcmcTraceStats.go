// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

func mcmcTraceStats(inFile string, outFile string, hdiProportion float64, burnIn int, parameterName string) {
	var t numbers.McmcTrace
	var err error
	t = numbers.ReadMcmcTrace(inFile, parameterName)
	numbers.DiscardBurnIn(t, burnIn)
	s, e := numbers.HighestDensityInterval(t, hdiProportion)
	m := numbers.MeanMcmcTrace(t)

	out := fileio.EasyCreate(outFile)
	defer out.Close()

	_, err = fmt.Fprintf(out, "Posterior mean: %f.\n", m)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "HDI Credible Interval. Proportion: %f. Start: %f. End: %f.\n", hdiProportion, s, e)
}

func usage() {
	fmt.Print(
		"mcmcTraceStats - Returns summary statistics on an MCMC trace file produced by selectionMCMC.\n" +
			"Usage:\n" +
			"mcmcTraceStats trace.txt out.txt\n")
	flag.PrintDefaults()
}


func main() {
	var expectedNumArgs int = 2
	var hdiProportion *float64 = flag.Float64("hdiProportion", 0.95, "Set the proportion of density contained in the hdi credible interval.")
	var burnIn *int = flag.Int("burnIn", 0, "Set the number of initial iterations to discard as burn-in from the trace.")
	var parameterName *string = flag.String("parameterName", "Mu", "Set the name of the parameter in the Mcmc trace to analyze.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	log.SetFlags(0)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	mcmcTraceStats(inFile, outFile, *hdiProportion, *burnIn, *parameterName)
}
