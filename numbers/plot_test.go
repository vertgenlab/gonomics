package numbers

import (
	"testing"
	"strings"
	"os"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
)



var PlotTests = []struct {
	f func(float64) float64
	left float64
	right float64
	bins int
	outFile string
	answerX []float64
	answerFofX []float64
	}{
	{func(x float64) float64 {return x}, 0.0, 3.0, 3, "test1.txt", []float64{0.0, 1.0, 2.0}, []float64{0.0, 1.0, 2.0}},
}

//TestPlot writes Plot results to a file, but then reads the file back in to confirm that they match expectations.
//Then it deletes the testing output temp file.
func TestPlot(t *testing.T) {
	var i int
	for _, test := range PlotTests {
		Plot(test.f, test.left, test.right, test.bins, test.outFile)
		file := fileio.EasyOpen(test.outFile)
		i = 0
		var line string
		var doneReading, firstTime bool = false, true
		for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
			if !firstTime {
					words := strings.Split(line, "\t")
				if common.StringToFloat64(words[0]) != test.answerX[i] {
					t.Errorf("TestPlot for X: %f returned: %f.", test.answerX[i], common.StringToFloat64(words[0]))
				}
				if common.StringToFloat64(words[1]) != test.answerFofX[i] {
					t.Errorf("TestPlot for f(X): %f returned: %f.", test.answerFofX[i], common.StringToFloat64(words[1]))
				}
				i++
			} else {
				firstTime = false
			}
		}
		file.Close()
		os.Remove(test.outFile)
	}
}