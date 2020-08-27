package simpleGraph

/*
import (
	"github.com/vertgenlab/gonomics/giraf"
	"sync"
	"testing"
)

func BenchmarkWriteToFileHandle(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var wg sync.WaitGroup
	for n := 0; n < b.N; n++ {
		reader := giraf.GoReadToChan("testdata/pairedTest.giraf")
		wg.Add(1)
		go giraf.GirafChanToFile("testdata/fileHandle.giraf", reader, &wg)
		wg.Wait()
	}
}

func BenchmarkStringBuilder(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var wg sync.WaitGroup
	for n := 0; n < b.N; n++ {
		reader := giraf.GoReadToChan("testdata/pairedTest.giraf")
		wg.Add(1)
		go WriteSimpleGiraf("testdata/stringBuilder.giraf", reader, &wg)
		wg.Wait()
	}
}
*/
