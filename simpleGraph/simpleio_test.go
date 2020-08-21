package simpleGraph

import (
	"testing"
	"os"
	"log"
	"runtime/pprof"
	"runtime"
	"flag"
)

var genomeGraphs = []struct {
	file string
}{
	{"testdata/rabsTHREEbepa_chain.gg"},
}
var cpuprofile string =  "simpleioProf/simpleio.ReaderCPU.prof"//flag.String("cpuprofile", "", "write cpu profile to `file`")
var memprofile string= "simpleioProf/simpleio.ReaderMEM.prof"//flag.String("memprofile", "", "write memory profile to `file`")
/*
func BenchmarkSimplePoolReader(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for _, test := range genomeGraphs {
		for n := 0; n < b.N; n++ {
			SimplyRead(test.file)
		}
		
	}	
}

func BenchmarkEasyReader(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for _, test := range  genomeGraphs {
		for n := 0; n < b.N; n++ {
			Read(test.file)
		}
	}	
}*/

func BenchmarkNewReader(b *testing.B) {
	flag.Parse()
    if cpuprofile != "" {
        f, err := os.Create(cpuprofile)
        if err != nil {
            log.Fatal("could not create CPU profile: ", err)
        }
        defer f.Close() // error handling omitted for example
        if err := pprof.StartCPUProfile(f); err != nil {
            log.Fatal("could not start CPU profile: ", err)
        }
        defer pprof.StopCPUProfile()
    }
	for _, test := range genomeGraphs {
		graph := SimplyRead(test.file)
		Write("testdata/check_rabsTHREEbepa.gg", &graph)
	}
	if memprofile != "" {
        f, err := os.Create(memprofile)
        if err != nil {
            log.Fatal("could not create memory profile: ", err)
        }
        defer f.Close() // error handling omitted for example
        runtime.GC() // get up-to-date statistics
        if err := pprof.WriteHeapProfile(f); err != nil {
            log.Fatal("could not write memory profile: ", err)
        }
    }

}