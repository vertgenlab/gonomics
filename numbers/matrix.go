package numbers

//Rref returns an input matrix in row reduced echelon form using Gaussian elimination.
func Rref(m [][]float64) [][]float64 {
	//first we make a copy of the input matrix to serve as our answer. This preserves the input matrix in its original form.
	mCopy := make([][]float64, len(m))
	for i := range m {
		mCopy[i] = make([]float64, len(m[i]))
		copy(mCopy[i], m[i])
	}

	var subtractFactor, entry float64
	var leadingCoefficient int = 0
	var i, j int
	for row := 0; row < len(mCopy); row++ {
		if leadingCoefficient >= len(mCopy[0]) {
			return mCopy
		}
		i = row
		for mCopy[i][leadingCoefficient] == 0 {
			i++
			if len(mCopy) == i {
				i = row
				leadingCoefficient++
				if leadingCoefficient == len(mCopy[0]) {
					return mCopy
				}
			}
		}
		mCopy[i], mCopy[row] = mCopy[row], mCopy[i]
		multFactor := 1 / mCopy[row][leadingCoefficient]
		for j = range mCopy[row] {
			mCopy[row][j] *= multFactor
		}
		for i = 0; i < len(mCopy); i++ {
			if i != row {
				subtractFactor = mCopy[i][leadingCoefficient]
				for j, entry = range mCopy[row] {
					mCopy[i][j] -= entry * subtractFactor
				}
			}
		}
		leadingCoefficient++
	}
	return mCopy
}