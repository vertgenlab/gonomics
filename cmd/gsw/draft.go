package main

/*
func GSWAlignerOne(ref *SimpleGraph, input string, output string) {
	var tileSize int = 12
	var stepSize int = 1

	outFile, _ := os.Create(output)

	defer outFile.Close()

	header := NodesHeader(ref.Nodes)


	sam.WriteHeaderToFileHandle(outFile, header)
	file := fileio.EasyOpen(input)
	defer file.Close()
	//runtime.GOMAXPROCS(runtime.NumCPU())

	threads := runtime.NumCPU()


	log.Printf("Indexing the genome...\n")
	ham5 := IndexGenomeDev(ref.Nodes, tileSize, stepSize)
	m, trace := SwMatrixSetup(10000)
	log.Printf("Aligning reads...\n")
	//queue := make(chan fastq.Fastq)

	tasks := make(chan *fastq.Fastq, 1000000)


	for workers := 0; workers < threads; workers++ {
		go func() {

			for {
				query, ok := <-tasks
				if !ok {
					//wg.Done()
					return
				}
				routinesGenomeGraph(ref, query, ham5, tileSize, m, trace)
				//mappedRead = GraphSmithWaterman(ref, query, ham5, tileSize, m, trace)
				//sam.WriteAlnToFileHandle(outFile, answer)
				//MapFastq(ref, query, seedLen, m, mSet, trace, outFile)
			}

		}()
	}

	var fq *fastq.Fastq
	var done bool
	//var i int

	//start implementation of goroutines
	var wg sync.WaitGroup
	inputChan := make(chan *fastq.Fastq, 10)
	results := make(chan *sam.SamAln, 10)
	var numReads int64 = 0
	//for x := 0; x < concurrency;x++ {
	//	go GoroutinesGenomeGraph(ref, inputChan, ham5, tileSize, m, trace, results)
	//}
	//allocate
	go func() {
		for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
			numReads++
			inputChan <- fq
		}
		close(inputChan)
	}()
	d := make(chan bool)

	go func(chan bool) {
		for j := range results {
			sam.WriteAlnToFileHandle(outFile, j)
		}
		d <- true
	}(d)

	var i int64
	//create pool
	go func() {
		for i =0; i < numReads; i++ {
			wg.Add(1)
			//go devGoroutinesGenomeGraph(ref, inputChan, ham5, tileSize, m, trace, results, &wg)
		}
		wg.Wait()
		close(results)
	}()


	<-d


	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		//mappedRead = GraphSmithWaterman(ref, fq, ham5, tileSize, m, trace)
		//sam.WriteAlnToFileHandle(outFile, mappedRead)

		//tasks <- fq

		go routinesGenomeGraph(ref, fq, ham5, tileSize, m, trace, outFile)

	}

	close(tasks)
	wg.Wait()
}

func GSWsBatch(ref *SimpleGraph, input string, output string, groupSize int) {
	//setting up dictionary/hash map look up, genome indexing
	var tileSize int = 15
	var stepSize int = 14
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeDev(ref.Nodes, tileSize, stepSize)
	m, trace := SwMatrixSetup(10000)

	//handling fastq file and creating sam file to write to output
	file := fileio.EasyOpen(input)
	defer file.Close()
	out, _ := os.Create(output)
	defer out.Close()
	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(out, header)

	//variables used to read and process fastq files
	var fq *fastq.Fastq
	var groups []*fastq.Fastq = make([]*fastq.Fastq, 0)
	var done bool
	c := make(chan *sam.SamAln, groupSize)

	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {

		groups = append(groups, fq)
		if len(groups) > groupSize {
			//send off group
			log.Printf("Aligning %d reads...\n", groupSize)
			routinesGenomeGraph(ref, groups, seedHash, tileSize, m, trace, c, out, groupSize)
			//zero out
			groups = groups[0:0]
		}
	}
	//last case check length groups
	if len(groups) != 0 {
		routinesGenomeGraph(ref, groups, seedHash, tileSize, m, trace, c, out, groupSize)
	}
}
/*
func NoLimitTestGSWsBatch(ref *SimpleGraph, input string, output string) {
	//setting up dictionary/hash map look up, genome indexing
	var tileSize int = 12
	var stepSize int = 1
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeDev(ref.Nodes, tileSize, stepSize)
	m, trace := SwMatrixSetup(10000)

	//handling fastq file and creating sam file to write to output
	file := fileio.EasyOpen(input)
	defer file.Close()
	out, _ := os.Create(output)
	defer out.Close()
	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(out, header)

	//variables used to read and process fastq files
	var fq *fastq.Fastq

	var done bool
	c := make(chan *sam.SamAln)
	var numReads int = 0
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		numReads++
		go wrap(ref, fq, seedHash, tileSize, m, trace, c)
	}
	for i := 0; i < numReads; i++ {
		sam.WriteAlnToFileHandle(out, <-c)
		//log.Printf("%s\n", sam.SamAlnToString(<-c))
	}

}
//no limit on number of goroutines
func noLimitGSW(ref *SimpleGraph, input string, seedHash [][]*SeedBed, seedLen int, m [][]int64, trace [][]rune) {
	//setting up dictionary/hash map look up, genome indexing




	//handling fastq file and creating sam file to write to output
	file := fileio.EasyOpen(input)
	defer file.Close()
	//out, _ := os.Create(output)
	//defer out.Close()
	//header := NodesHeader(ref.Nodes)
	//sam.WriteHeaderToFileHandle(out, header)

	//variables used to read and process fastq files
	var fq *fastq.Fastq

	var done bool
	c := make(chan *sam.SamAln, 1000)
	log.Printf("Aligning reads...\n")
	var numReads int = 0
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		numReads++
		go wrap(ref, fq, seedHash, seedLen, m, trace, c)
	}
	for i := 0; i < numReads; i++ {
		//sam.WriteAlnToFileHandle(out, <-c)
		log.Printf("%s\n", sam.SamAlnToString(<-c))
	}

}

//func GSWsBatch(ref *SimpleGraph, input string, divide int) {
func devGSWsBatch(ref *SimpleGraph, input string, seedHash [][]*SeedBed, seedLen int, m [][]int64, trace [][]rune, groupSize int) {
	c := make(chan *sam.SamAln, groupSize)
	//var tileSize int = 12
	//var stepSize int = 1
	//log.Printf("Indexing the genome...\n")
	//seedHash := IndexGenomeDev(ref.Nodes, tileSize, stepSize)
	//m, trace := SwMatrixSetup(10000)


	var fq *fastq.Fastq

	var groups []*fastq.Fastq = make([]*fastq.Fastq, 0)
	var done bool
	//fastq input
	file := fileio.EasyOpen(input)
	defer file.Close()

	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {

		groups = append(groups, fq)
		if len(groups) > groupSize {
			//send off group
			log.Printf("Aligning %d reads...\n", groupSize)
			devGoroutinesGenomeGraph(ref, groups, seedHash, seedLen, m, trace, c)
			//zero out
			groups = groups[0:0]
		}
	}
	//last case check length groups
	if len(groups) != 0 {
		devGoroutinesGenomeGraph(ref, groups, seedHash, seedLen, m, trace, c)
	}
}*/
