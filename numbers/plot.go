package numbers

import (
	"github.com/vertgenlab/gonomics/fileio"
	"fmt"
)

//Plot returns the values of a function(float64) float64 between a left and right bound for a specified number of bins.
//the answer is written to a file as a CSV for subsequent visualization.
//Half-closed interval [left, right)
func Plot(f func(float64) float64, left float64, right float64, bins int, outFile string) {
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	var step float64 = (right - left) / float64(bins)
	var current float64 = left

	fmt.Fprintf(out, "X\tf(X)\n")

	for i := 0; i < bins; i++ {
		fmt.Fprintf(out, "%f\t%f\n", current, f(current))
		current = current + step
	}
}