package numbers

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"math"
	"sort"
	"strings"
)

//McmcTrace is a general struct for Mcmc trace output. Used for discarding burn-in and calculating the mean and credible interval.
type McmcTrace struct {
	Parameter []float64 //Parameter state, where Parameter[i] is the value of Parameter in the ith iteration.
}

func ReadMcmcTrace(inFile string, parameterName string) McmcTrace {
	var curr string
	var words []string
	var doneReading bool
	var currParam float64
	var t McmcTrace = McmcTrace{make([]float64, 0)}

	in := fileio.EasyOpen(inFile)
	ParameterIndex := parseMcmcTraceHeader(in, parameterName)

	for curr, doneReading = fileio.EasyNextLine(in); !doneReading; curr, doneReading = fileio.EasyNextLine(in) {
		words = strings.Split(curr, "\t")
		currParam = common.StringToFloat64(words[ParameterIndex])
		t.Parameter = append(t.Parameter, currParam)
	}

	return t
}

//parseMcmcTraceHeader is a helper function of ReadMcmcTrace that finds the column of the parameter to be analyzed.
func parseMcmcTraceHeader(in *fileio.EasyReader, parameterName string) int {
	var headerLine string
	var doneReading bool
	var words []string
	var ParameterIndex int = -1 //stores the column index that stores the parameter to analyze.

	headerLine, doneReading = fileio.EasyNextLine(in)
	if doneReading {
		log.Fatalf("Empty trace file.")
	}
	words = strings.Split(headerLine, "\t")
	if words[0] != "Iteration" {
		log.Fatalf("Improperly formatted MCMC trace file. Trace files generated with selectionMCMC will begin with a header line starting with 'Iteration'.")
	}
	for i := 1; i < len(words); i++ {
		if words[i] == parameterName {
			ParameterIndex = i
		}
	}
	if ParameterIndex == -1 {
		log.Fatalf("No column with the input parameterName, %s, is found in the trace file.\n", parameterName)
	}
	return ParameterIndex
}

//DiscardBurnIn will remove the the first i values in an McmcTrace, where i is equal to the input value burnIn.
func DiscardBurnIn(t McmcTrace, burnIn int) {
	t.Parameter = t.Parameter[burnIn:]
}

//HighestDensityInterval returns the HDI credible interval for an input McmcTrace struct. Proportion is the proportion of iterations in the credible interval (ex. 0.95 for a 95% credible interval).
func HighestDensityInterval(t McmcTrace, proportion float64) (float64, float64) {
	var minStart, minEnd float64
	var tmp []float64 = make([]float64, len(t.Parameter))
	copy(tmp, t.Parameter)
	sort.Float64s(tmp) //sorts low to high

	var pIndex int = int(math.Ceil(proportion*float64(len(tmp)))) - 1 //returns the index of the highest value at the left-tailed credible interval.
	minStart, minEnd = tmp[0], tmp[pIndex]
	var minDistance float64 = minEnd - minStart
	var currDistance float64

	for i := 1; i < len(tmp)-pIndex; i++ { //loop through all possible 95% credible intervals
		currDistance = tmp[pIndex+i] - tmp[i]
		if currDistance < minDistance {
			minStart, minEnd = tmp[i], tmp[pIndex+i]
			minDistance = currDistance
		}
	}

	return minStart, minEnd
}

//MeanMcmcTrace returns the mean value of the posterior distribution estimated by an McmcTrace.
func MeanMcmcTrace(t McmcTrace) float64 {
	return AverageFloat64(t.Parameter)
}
