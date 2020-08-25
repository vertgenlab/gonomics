package giraf

import (
	"testing"
	"sync"
)

func TestReadNewWrite(t *testing.T) {
	Read("testdata/stringBuilder.giraf")
} 

func BenchmarkWriteToFileHandle(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var wg sync.WaitGroup
	for n := 0; n < b.N; n++ {
		reader := GoReadToChan("testdata/pairedTest.giraf")
		wg.Add(1)
		go GirafChanToFile("testdata/fileHandle.giraf", reader, &wg)
		wg.Wait()
	}
}

func BenchmarkStringBuilder(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var wg sync.WaitGroup
	for n := 0; n < b.N; n++ {
		reader := GoReadToChan("testdata/pairedTest.giraf")
		wg.Add(1)
		go WriteSimpleGiraf("testdata/stringBuilder.giraf", reader, &wg)
		wg.Wait()
	}
}
