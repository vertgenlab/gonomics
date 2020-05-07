package main
import(

	"log"
	"sync"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"path/filepath"
	"time"
	"strings"
)

func GswToGirafPair(ref *simpleGraph.SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64) {
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := simpleGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup

	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	girafPipe := make(chan *giraf.GirafPair, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)
	log.Printf("Scoring matrix used:\n%s\n", simpleGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)

	log.Printf("Aligning %s and %s to genome graph...", strings.Split(filepath.Base(readOne), ".")[0], strings.Split(filepath.Base(readTwo), ".")[0])
	start := time.Now()
	for i := 0; i < threads; i++ {
		go simpleGraph.RoutineFqPairToGiraf(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, girafPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go giraf.GirafPairChanToFile(output, girafPipe, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(girafPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}


func WrapGirafLiftoverToSam(ref *simpleGraph.SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64, header *sam.SamHeader) {
	log.SetFlags(log.Ldate | log.Ltime)
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := simpleGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	//log.Printf("Setting up read and write channels...\n\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)

	log.Printf("Scoring matrix used:\n%s\n", simpleGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning sequence to genome graph...")
	start := time.Now()
	for i := 0; i < threads; i++ {
		go simpleGraph.RoutineGirafToSam(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}