package wig

// Smooth performs moving average smoothing on an input Wig struct using a user-specified windowSize.
func Smooth(w Wig, windowSize int, missing float64) Wig {
	var j, k, midPoint int
	var sum float64
	var missingFound bool = false
	//first we make a fresh copy of the input wig to serve as our output. Default value is "missing".
	var answer = Wig{StepType: "fixedStep", Chrom: w.Chrom, Start: w.Start, Step: w.Step, Span: w.Span, Values: make([]float64, len(w.Values))}
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

// SmoothSlice performs moving average smoothing on an input slice of Wig structs using a user-specified windowSize.
func SmoothSlice(w []Wig, windowSize int, missing float64) []Wig {
	var answer = make([]Wig, len(w))
	for i := range w {
		answer[i] = Smooth(w[i], windowSize, missing)
	}
	return answer
}
