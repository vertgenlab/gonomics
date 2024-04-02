package matrix

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math"
	"strings"
)

type Matrix [][]float64

func Print(m Matrix) {
	for i := range m {
		fmt.Println(m[i])
	}
}

func Minimum(m Matrix) float64 {
	var ans float64 = m[0][0]
	var j int
	for i := range m {
		for j = range m[i] {
			if m[i][j] < ans {
				ans = m[i][j]
			}
		}
	}
	return ans
}

func AppendRow(m *Matrix, newRow []float64) {
	nRow, nCol := Shape(*m)
	if nCol != len(newRow) && nRow != 0 {
		log.Fatalf("New row doesn't have the same number of elements as the matrix.")
	}
	*m = append(*m, newRow)
}

func SubsetRows(m Matrix, idx []int) Matrix {
	n := Matrix{}
	for _, i := range idx {
		n = append(n, m[i])
	}
	return n
}

// Where returns a [][]int where i,j is an the index of the location in the matrix that satisfies the input expression
func Where(m Matrix, expression string, value float64) [][]int {
	var j int
	ans := make([][]int, 2)
	switch expression {
	case "=":
		for i := range m {
			for j = range m[i] {
				if m[i][j] == value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	case ">":
		for i := range m {
			for j = range m[i] {
				if m[i][j] > value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	case ">=":
		for i := range m {
			for j = range m[i] {
				if m[i][j] >= value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	case "<":
		for i := range m {
			for j = range m[i] {
				if m[i][j] < value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	case "<=":
		for i := range m {
			for j = range m[i] {
				if m[i][j] <= value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	case "!=":
		for i := range m {
			for j = range m[i] {
				if m[i][j] != value {
					ans[0] = append(ans[0], i)
					ans[1] = append(ans[1], j)
				}
			}
		}
	}
	return ans
}

func Maximum(m Matrix) float64 {
	var ans float64 = m[0][0]
	var j int
	for i := range m {
		for j = range m[i] {
			if m[i][j] > ans {
				ans = m[i][j]
			}
		}
	}
	return ans
}

// AllSums takes a matrix returns a matrix where the first row is the sums of all the rows and second row is the sums of all the columns. All ordered by index.
func AllSums(m Matrix) Matrix {
	_, c := Shape(m)
	n := Matrix{}
	n = append(n, []float64{})
	for i := range m {
		n[0] = append(n[0], SumRow(m, i))
	}
	n = append(n, []float64{})
	for i := 0; i < c; i++ {
		n[1] = append(n[1], SumColumn(m, i))
	}
	return n
}

func Average(m Matrix) float64 {
	var total float64
	r, c := Shape(m)
	for i := range m {
		total += SumRow(m, i)
	}
	return total / float64(r*c)
}

func SameShape(a, b Matrix) bool {
	rowA, colA := Shape(a)
	rowB, colB := Shape(b)
	if rowA == rowB && colA == colB {
		return true
	} else {
		return false
	}
}

func AppendColumn(m Matrix, values []float64) {
	nRow, _ := Shape(m)
	if nRow != len(values) {
		log.Fatalf("New column doesn't have the same number of rows as the matrix")
	}
	for i := range values {
		m[i] = append(m[i], values[i])
	}
}

func Combine(a, b Matrix, operation rune) Matrix {
	var j int
	if !SameShape(a, b) {
		log.Fatalf("Matricies are not of the same size")
	}
	n := Initilize(Shape(a))
	switch operation {
	case '+':
		for i := range a {
			for j = range a[i] {
				n[i][j] = a[i][j] + b[i][j]
			}
		}
	case '*':
		for i := range a {
			for j = range a[i] {
				n[i][j] = a[i][j] * b[i][j]
			}
		}
	case '/':
		for i := range a {
			for j = range a[i] {
				n[i][j] = a[i][j] / b[i][j]
			}
		}
	case '-':
		for i := range a {
			for j = range a[i] {
				n[i][j] = a[i][j] - b[i][j]
			}
		}
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

func SumRow(mat Matrix, idx int) float64 {
	var ans float64
	for i := range mat[idx] {
		ans += mat[idx][i]
	}
	return ans
}

func Diags(f []float64) Matrix {
	m := Initilize(len(f), len(f))
	for i := range f {
		m[i][i] = f[i]
	}
	return m
}

func SumColumn(mat Matrix, idx int) float64 {
	var ans float64
	for row := range mat {
		ans += mat[row][idx]
	}
	return ans
}

func Sum(m Matrix) float64 {
	var ans float64
	for i := range m {
		ans += SumRow(m, i)
	}
	return ans
}

func Flatten(m Matrix) Matrix {
	var n Matrix
	n = append(n, []float64{})
	var j int
	for i := range m {
		for j = range m[i] {
			n[0] = append(n[0], m[i][j])
		}
	}
	return n
}

func GetColumn(mat Matrix, idx int) []float64 {
	var col []float64
	for r := range mat {
		col = append(col, mat[r][idx])
	}
	return col
}

// Shape returns rows x columns
func Shape(m Matrix) (int, int) {
	return len(m), len(m[0])
}

func Transpose(m Matrix) Matrix {
	r, c := Shape(m)
	n := Initilize(c, r)
	for i := range m[0] {
		n[i] = GetColumn(m, i)
	}
	return n
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

func Read(in, delim string) Matrix {
	var m Matrix
	var j int
	inMat := fileio.Read(in)
	for i := range inMat {
		vals := strings.Split(inMat[i], delim)
		m = append(m, []float64{})
		for j = range vals {
			m[i] = append(m[i], parse.StringToFloat64(vals[j]))
		}
	}
	return m
}

func rowToString(r []float64, delim string) string {
	var slc []string
	for i := range r {
		slc = append(slc, fmt.Sprintf("%f", r[i]))
	}
	return strings.Join(slc, delim)
}

func Write(file *fileio.EasyWriter, m Matrix, delim string) {
	for i := range m {
		fileio.WriteToFileHandle(file, rowToString(m[i], delim))
	}
}

func IsDoubleStochastic(m Matrix) bool {
	r, c := Shape(m)
	if r != c {
		log.Fatalf("matrix isn't symmetric. Rows: %d Columns: %d\n", r, c)
	}
	for i := range m {
		if math.Abs(1-SumRow(m, i)) > 1e-6 {
			fmt.Printf("Row %d sum: %f\n", i, SumRow(m, i))
			return false
		}
		if math.Abs(1-SumColumn(m, i)) > 1e-6 {
			fmt.Printf("Column %d sum: %f\n", i, SumColumn(m, i))
			return false
		}
	}
	return true
}
