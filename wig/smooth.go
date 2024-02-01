package wig

// Smooth performs moving average smoothing on an input Wig struct using a user-specified windowSize and a value that represents missing data in the Wig.
func Smooth(w Wig, windowSize int, missing float64) Wig {
	var j, k, midPoint int
	var sum float64
	var missingFound bool = false
	//first we make a fresh copy of the input wig to serve as our output. Default value is "missing".
	var answer = Wig{StepType: "fixedStep", Chrom: w.Chrom, Start: w.Start, Step: w.Step, Span: w.Span, DefaultValue: w.DefaultValue, Values: make([]float64, len(w.Values))}
	copy(answer.Values, w.Values)
	for i := range answer.Values {
		answer.Values[i] = missing
	}

	for j = 0; j < len(w.Values)-windowSize; j++ {
		missingFound = false
		sum = 0.0
		for k = j; k < j+windowSize; k++ {
			if w.Values[k] == missing {
				missingFound = true
				break
			}
			sum += w.Values[k]
		}
		if !missingFound {
			midPoint = (j + j + windowSize) / 2
			answer.Values[midPoint] = sum / float64(windowSize)
		}
	}
	return answer
}

// SmoothMap performs moving average smoothing on an input map of Wig structs using a user-specified windowSize and a value that represents missing data in the Wig.
func SmoothMap(w map[string]Wig, windowSize int, missing float64) map[string]Wig {
	var answer = make(map[string]Wig)
	for currKey := range w {
		answer[currKey] = Smooth(w[currKey], windowSize, missing)
	}
	return answer
}
