package wig

// Smooth performs moving average smoothing on an input Wig struct using a user-specified windowSize.
func Smooth(w Wig, windowSize int) Wig {
	var j, k, midPoint int
	var sum float64
	//first we make a fresh copy of the input wig to serve as our output.
	var answer = Wig{StepType: "fixedStep", Chrom: w.Chrom, Start: w.Start, Step: w.Step, Span: w.Span, Values: make([]float64, len(w.Values))}
	copy(answer.Values, w.Values)
	for j = 0; j < len(w.Values)-windowSize; j++ {
		sum = 0.0
		for k = j; k < j+windowSize; k++ {
			sum += w.Values[k]
		}
		midPoint = (j + j + windowSize) / 2
		answer.Values[midPoint] = sum / float64(windowSize)
	}
	return answer
}

// SmoothSlice performs moving average smoothing on an input slice of Wig structs using a user-specified windowSize.
func SmoothSlice(w []Wig, windowSize int) []Wig {
	var answer = make([]Wig, len(w))
	for i := range w {
		answer[i] = Smooth(w[i], windowSize)
	}
	return answer
}
