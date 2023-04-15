package fileio

import (
	"testing"
)

func BenchmarkGunzipReaderReg(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewGunzipReader("testdata/big.fa")
		for _, done := GunzipLine(reader); !done; _, done = GunzipLine(reader) {
			// Nothing to assign, testing pure reading of the file
		}
	}
}

func BenchmarkGunzipReaderGz(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		reader := NewGunzipReader("testdata/big.fa.gz")
		for _, done := GunzipLine(reader); !done; _, done = GunzipLine(reader) {
			// Nothing to assign, testing pure reading of the file
		}
	}
}
