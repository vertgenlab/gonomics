package giraf

import (
	"testing"
	//"github.com/vertgenlab/gonomics/dna"
)

var warrior *Giraf = Read("testdata/fileHandle.giraf")[0]

func BenchmarkStringBuilder(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		//for _, warrior := range gsw {
		ToString(warrior)
		//}

	}
}

func BenchmarkGirfToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		//	for _, warrior := range gsw {
		GirafToString(warrior)

	}
}
