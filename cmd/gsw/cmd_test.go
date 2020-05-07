package main
/*
import(
	"testing"
)

func TestExcute(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 0
	var cpus int = 4
	var scoreMatrix = HumanChimpTwoScoreMatrix
	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/bigGenome.sg")
	chrSize := chromInfo.ReadToSlice("testdata/bigGenome.sizes")
	header := sam.ChromInfoSamHeader(chrSize)
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	readOne := "testdata/simReads_R1.fastq"
	readTwo := "testdata/simReads_R2.fastq"
	os.Remove(readOne)
	os.Remove(readTwo)
	fastq.WritePair(readOne, readTwo, simReads)
	GSWsPair(genome, readOne, readTwo, "testdata/test.sam", cpus, tileSize, stepSize, scoreMatrix, header)
	samfile, _ := sam.Read("testdata/test.sam")
	for _, samline := range samfile.Aln {
		log.Printf("%s\n", ViewGraphAlignment(samline, genome))
	}
	os.Remove("testdata/test.sam")
}

func TestPabBioGraph(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix

	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/bigGenome.sg")

	//genome, _ := Read("testdata/rabsBepaChrI.gg")

	log.Printf("Indexing the genome...\n")
	tiles := indexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.FastqBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)

	go fastq.ReadBigToChan("testdata/simReads.fq", fastqPipe)
	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "testdata/test.sam", nil, &writerWaiter)

	log.Printf("Starting alignment worker...\n")
	time.Sleep(5 * time.Second)

	f, err := os.Create("testdata/cpuprofile.data")
	common.ExitIfError(err)
	defer f.Close()
	err = pprof.StartCPUProfile(f)
	common.ExitIfError(err)

	start := time.Now()
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswWorkerMemPool(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, samPipe, &workerWaiter)
	}

	workerWaiter.Wait()
	stop := time.Now()
	pprof.StopCPUProfile()
	close(samPipe)

	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	duration := stop.Sub(start)

	samfile, _ := sam.Read("testdata/test.sam")
	for _, samline := range samfile.Aln {
		log.Printf("%s\n", ViewGraphAlignment(samline, genome))
	}
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
	//strictCheckFile("testPabBioGraph.sam", genome)

	//log.Printf("Passed alignment check!!!\n")
}

func TestGirafLiftoverToSam(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1000
	var readLength int = 150
	var mutations int = 0
	var cpus int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/bigGenome.sg")
	chrSize := chromInfo.ReadToSlice("testdata/bigGenome.sizes")
	header := sam.ChromInfoSamHeader(chrSize)
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	readOne := "testdata/simReads_R1.fastq"
	readTwo := "testdata/simReads_R2.fastq"
	os.Remove(readOne)
	os.Remove(readTwo)
	fastq.WritePair(readOne, readTwo, simReads)
	WrapGirafLiftoverToSam(genome, readOne, readTwo, "testdata/liftoverToSam.sam", cpus, tileSize, stepSize, scoreMatrix, header)
}

func TestExcuteGiraf(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 0
	var cpus int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/bigGenome.sg")

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	readOne := "testdata/simReads_R1.fastq"
	readTwo := "testdata/simReads_R2.fastq"
	os.Remove(readOne)
	os.Remove(readTwo)
	fastq.WritePair(readOne, readTwo, simReads)
	GswToGirafPair(genome, readOne, readTwo, "/dev/stdout", cpus, tileSize, stepSize, scoreMatrix)
}

func TestGirafGSW(t *testing.T) {
	var output string = "testdata/giraf.tsv"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	//genome := Read("testdata/bepaRabsSV.gg")
	genome := Read("testdata/bigGenome.sg")
	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)

	log.Printf("Making sam channel...\n")
	girafPipe := make(chan *giraf.GirafPair, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	tiles := simpleGraph.IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	go fastq.ReadPairBigToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go simpleGraph.RoutineFqPairToGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}
	go giraf.GirafPairChanToFile(output, girafPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	//samfile, _ := sam.Read(output)
	//for _, samline := range samfile.Aln {
	///	log.Printf("%s\n", ViewGraphAlignment(samline, genome))
	//}
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
	//os.Remove(output)
}*/