package popgen

import "testing"

func BenchmarkAfsF1e7(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-7)
	}
}

func BenchmarkAfsF1e6(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-6)
	}
}

func BenchmarkAfsF1e5(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-5)
	}
}

func BenchmarkAfsF1e4(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-4)
	}
}
