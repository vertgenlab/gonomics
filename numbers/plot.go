package numbers

import (
	"fmt"

	"github.com/vertgenlab/gonomics/fileio"
)

//Plot returns the values of a function(float64) float64 between a left and right bound for a specified number of bins.
//the answer is written to a file as a CSV for subsequent visualization.
//Half-closed interval [left, right)
func Plot(f func(float64) float64, left float64, right float64, bins int, outFile string) {
	out := fileio.EasyCreate(outFile)
	var step float64 = (right - left) / float64(bins)
	var current float64 = left

	fmt.Fprintf(out, "X\tf(X)\n")

	for i := 0; i < bins; i++ {
		fmt.Fprintf(out, "%f\t%f\n", current, f(current))
		current = current + step
	}
	out.Close()
}

//PlotBinomCoefficient writes binomial coefficients (n choose k) from k=1 to k=n-1 to an output file for downstream visualization.
func PlotBinomCoefficient(n int, outFile string) {
	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "i\tProbability\n")

	for i := 1; i < n; i++ { //all possible choose, not including 0 and n as i arguments.
		fmt.Fprintf(out, "%v\t%v\n", i, BinomCoefficientLog(n, i))
	}
	out.Close()
}
