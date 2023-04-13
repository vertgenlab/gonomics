package common

import (
	"testing"

	"golang.org/x/exp/constraints"
	"golang.org/x/exp/slices"
)

func intContains(s []int, v int) bool {
	for i := range s {
		if s[i] == v {
			return true
		}
	}
	return false
}

func genericContains[E comparable](s []E, v E) bool {
	for i := range s {
		if s[i] == v {
			return true
		}
	}
	return false
}

func BenchmarkIntContains(b *testing.B) {
	s := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	for i := 0; i < b.N; i++ {
		intContains(s, 9)
	}
}

func BenchmarkGenericContains(b *testing.B) {
	s := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	for i := 0; i < b.N; i++ {
		genericContains(s, 9)
	}
}

func BenchmarkBuiltinContains(b *testing.B) {
	s := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	for i := 0; i < b.N; i++ {
		slices.Contains(s, 9)
	}
}

func intMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func genericMin[E constraints.Ordered](a, b E) E {
	if a < b {
		return a
	}
	return b
}

const (
	one int = 1
	two int = 2
)

func BenchmarkIntMin(b *testing.B) {
	for i := 0; i < b.N; i++ {
		intMin(one, two)
	}
}

func BenchmarkGenericMin(b *testing.B) {
	for i := 0; i < b.N; i++ {
		genericMin(one, two)
	}
}
