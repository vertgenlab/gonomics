package variant

import "testing"

const mag = 10

func BenchmarkAppend(b *testing.B) {
	ten := make([]int, 1 *mag)
	eighty := make([]int, 8 *mag)
	for i := 0; i < b.N; i++ {
		a := make([]int, 0, 10 *mag)
		a = append(a, ten...)
		a = append(a, eighty...)
		a = append(a, ten...)
	}
}

func BenchmarkCopy(b *testing.B) {
	ten := make([]int, 1 *mag)
	eighty := make([]int, 8 *mag)
	for i := 0; i < b.N; i++ {
		a := make([]int, 10 *mag)
		copy(a[:1*mag], ten)
		copy(a[1*mag:9*mag], eighty)
		copy(a[9*mag:], ten)
	}
}
