package numbers

//BinomialDistLog returns log(BinomialDist), where log is the natural logarithm.
//This is ideal for very small probabilities to avoid underflow. 
func BinomialDistLog(n int, k int, p float64) float64 {
	coefficient := BinomCoefficientLog(n, k)
	s := LogPowInt(p, k)
	f := LogPowInt(1.0-p, n-k)
	expression := MultiplyLog(s, f)
	return MultiplyLog(coefficient, expression)
}
