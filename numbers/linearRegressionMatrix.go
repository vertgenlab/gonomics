package numbers

import (
	"fmt"
	"log"
)

type Matrix [][]float64

func LinearRegression(m Matrix) (slope, itercept float64) {
	var x, y Matrix
	tempX := GetColumn(m, 0)
	for i := range tempX {
		x = append(x, []float64{1, tempX[i]})
	}
	tempY := GetColumn(m, 1)
	for i := range tempY {
		y = append(y, []float64{tempY[i]})
	}
	xT := Transpose(x)
	xTdotT := DotProduct(xT, x)
	xTdotTinverse := Inverse2x2Matrix(xTdotT)
	dotProductWithInverse := DotProduct(xTdotTinverse, xT)
	finalMatrix := DotProduct(dotProductWithInverse, y)
	return finalMatrix[1][0], finalMatrix[0][0]
}

func Transpose(m Matrix) Matrix {
	r, c := Shape(m)
	n := Initilize(c, r)
	for i := range m[0] {
		n[i] = GetColumn(m, i)
	}
	return n
}

func Initilize(rows, columns int) Matrix {
	mat := make([][]float64, rows)
	for r := range mat {
		mat[r] = make([]float64, columns)
	}
	return mat
}

// Shape returns rows x columns
func Shape(m Matrix) (int, int) {
	return len(m), len(m[0])
}

func GetColumn(mat Matrix, idx int) []float64 {
	var col []float64
	for r := range mat {
		col = append(col, mat[r][idx])
	}
	return col
}

func Inverse2x2Matrix(m Matrix) Matrix {
	n := Initilize(2, 2)
	n[0][0] = m[1][1]
	n[0][1] = -m[0][1]
	n[1][0] = -m[1][0]
	n[1][1] = m[0][0]
	return Scale(n, 1/(m[0][0]*m[1][1]-m[0][1]*m[1][0]), '*')
}

func Print(m Matrix) {
	for i := range m {
		fmt.Println(m[i])
	}
}

func DotProduct(a, b Matrix) Matrix {
	var i, j, k int
	var total float64
	rowA, colA := Shape(a)
	rowB, colB := Shape(b)

	if colA != rowB {
		log.Fatalf("Matrix multiplication not compatible for matricies of this size.")
	}

	result := Initilize(rowA, colB)

	for i = 0; i < rowA; i++ {
		for j = 0; j < colB; j++ {
			for k = 0; k < rowB; k++ {
				total += a[i][k] * b[k][j]
			}
			result[i][j] = total
			total = 0
		}
	}

	return result
}

func Scale(m Matrix, val float64, operation rune) Matrix {
	n := Initilize(Shape(m))
	switch operation {
	case '-':
		for i := range m {
			for j := range m[i] {
				n[i][j] = m[i][j] - val
			}
		}
	case '+':
		for i := range m {
			for j := range m[i] {
				n[i][j] = m[i][j] + val
			}
		}
	case '*':
		for i := range m {
			for j := range m[i] {
				n[i][j] = m[i][j] * val
			}
		}
	case '/':
		for i := range m {
			for j := range m[i] {
				n[i][j] = m[i][j] / val
			}
		}
	}
	return n
}
