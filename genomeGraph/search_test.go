package genomeGraph

import (
	"flag"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sync"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestRightDynamicAln(t *testing.T) {
	// Setup
	alpha := dna.StringToBases("AGTGGCTATCGTC") // Example sequence.
	beta := dna.StringToBases("GGCTAT")

	memory := NewMatrixPool()

	// Call DynamicAln function
	score, cig, endAlpha, endBeta := RightDynamicAln(alpha, beta, align.DefaultScoreMatrix, memory, -100, dynamicScoreKeeper{})
	settings := SearchSettings{
		scoreMatrix: align.DefaultScoreMatrix,
		gapPenalty:  -100,
		tileSize:    32,
		extension:   150,
	}
	pool := NewCache(len(alpha) + len(beta))
	visited := pool.Get().(*SearchCache)
	newScore, newCig, newTEnd, newQEnd := RightGraphAlign(alpha, beta, *visited, settings)
	// Expected results
	// var expectedScore int64 = 273                                                                       // Example expected score.
	// var expectedCigar []cigar.ByteCigar = []cigar.ByteCigar{{RunLen: 6, Op: 'M'}, {RunLen: 3, Op: 'D'}} // Example expected CIGAR.
	// var expectedEndAlpha, expectedEndBeta int = 9, 6

	// Log results
	t.Logf("Alignment Score: %d\nAlignment: %s\n", newScore, cigar.ByteCigarToString(newCig))

	// Check results
	if score != newScore {
		t.Errorf("Error: Score %v != %v\n", score, newScore)
	}
	if cigar.ByteCigarToString(cig) != cigar.ByteCigarToString(newCig) {
		t.Errorf("Error: CIGAR %+v != %+v\n", cig, newCig)
	}
	if endAlpha != newTEnd || endBeta != newQEnd {
		t.Errorf("Error: Ending indices (%d, %d) != (%d, %d) do not match expected values\n", newTEnd, newQEnd, endAlpha, endBeta)
	}
}

func TestAlignGraphTraversal(t *testing.T) {
	// Setup a simple genome graph for testing
	genome := EmptyGraph()

	var n0, n1, n2 Node
	var e0, e1 Edge
	/*
		Graph nodes: 		AAAAAT -> ATCG -> CGTTTTTT
		Input Read: 		ATATCGCGT
		Output Expected:	AT->ATCG->CGT
	*/
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("AAAAAT"),
	}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("ATCG"),
	}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("CGTTTTTT"),
	}

	// Make Edges
	e0 = Edge{
		Dest: &n1,
		Prob: 1,
	}

	e1 = Edge{
		Dest: &n2,
		Prob: 1,
	}
	// Define Paths
	n0.Next = append(n0.Next, e0, e1)
	n1.Prev = append(n1.Prev, Edge{Dest: &n0, Prob: 1})
	n2.Prev = append(n2.Prev, Edge{Dest: &n0, Prob: 1})

	genome.Nodes = append(genome.Nodes, n0, n1, n2)

	pool := &sync.Pool{
		New: func() interface{} {
			return &dnaPool{}
		},
	}

	memory := NewMatrixPool()

	read := dna.StringToBases("ATATCGCGT")

	dynamicScore := dynamicScoreKeeper{}
	scoreKeeper := scoreKeeper{}

	// expectedRightPath :=  []uint32{1,2}
	// Expected results for RightDirection
	expectedRightScore := int64(850) // Simplified expected score
	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 9, Op: 'M'}}
	expectedRightTargetEnd := 9
	expectedRightQueryEnd := 9

	rightAlign, rightScore, rightTargetEnd, rightQueryEnd, _ := RightAlignTraversal(&genome.Nodes[0], read, 0, []uint32{}, len(read), read, align.HumanChimpTwoScoreMatrix, memory, &scoreKeeper, &dynamicScore, pool)

	leftAlign, leftScore, leftTargetStart, leftQueryStart, _ := LeftAlignTraversal(&genome.Nodes[0], read, 0, []uint32{}, len(read), read, align.HumanChimpTwoScoreMatrix, memory, scoreKeeper, dynamicScore, pool)

	if rightScore != expectedRightScore {
		t.Errorf("RightDirection: expected score %d, got %d", expectedRightScore, rightScore)
	}
	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(expectedRightAlignment) {
		t.Errorf("RightDirection: expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
	}
	if expectedRightTargetEnd != rightTargetEnd || rightQueryEnd != expectedRightQueryEnd {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetEnd, rightQueryEnd, rightTargetEnd, rightQueryEnd)
	}
	settings := SearchSettings{
		scoreMatrix: align.HumanChimpTwoScoreMatrix,
		gapPenalty:  -100,
		tileSize:    32,
		extension:   len(read),
	}
	visited := NewCache(defaultMatrixSize)

	newRightScore, newRightAlign, newRightTargetEnd, newRightQueryEnd, _ := BreathSearchRight(&genome.Nodes[0], read, 0, []uint32{}, &visited, settings)
	newLeftScore, newLeftAlign, newLeftTargetStart, newLeftQueryStart, _ := BreathSearchLeft(&genome.Nodes[0], read, 0, []uint32{}, &visited, settings)
	if leftScore != newLeftScore {
		t.Errorf("LeftDirection: expected score %d, got %d", leftScore, newLeftScore)
	}
	if cigar.ByteCigarToString(leftAlign) != cigar.ByteCigarToString(newLeftAlign) {
		t.Errorf("LeftDirection: expected alignment %+v, got %+v", newRightAlign, newLeftAlign)
	}
	if rightScore != newRightScore {
		t.Errorf("RightDirection: expected score %d, got %d", newRightScore, rightScore)
	}
	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(newRightAlign) {
		t.Errorf("RightDirection: expected alignment %+v, got %+v", newRightAlign, rightAlign)
	}
	if newRightTargetEnd != rightTargetEnd || rightQueryEnd != newRightQueryEnd {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", newRightTargetEnd, newRightQueryEnd, rightTargetEnd, rightQueryEnd)
	}
	if leftTargetStart != newLeftTargetStart || leftQueryStart != newLeftQueryStart {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", newRightTargetEnd, newRightQueryEnd, rightTargetEnd, rightQueryEnd)
	}
}

var cpuprofile = flag.String("cpuprofile", "cpuprofile", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "memprofile", "write memory profile to `file`")

func BenchmarkGirafAlignment(b *testing.B) {
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	// var output string = "testdata/rabs_test.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1
	var readLength int = 150
	var mutations int = 5
	//var workerWaiter, writerWaiter sync.WaitGroup
	//var numWorkers int = 6
	//var scoreMatrix = align.HumanChimpTwoScoreMatrix

	genome := Read("testdata/bigGenome.sg")
	fqOne, fqTwo := "testdata/simReads_R1.fq", "testdata/simReads_R2.fq"

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)

	fastq.WritePair(fqOne, fqTwo, simReads)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	// go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	fqs := []fastq.FastqBig{}
	file := fileio.NewByteReader(fqOne)
	for read, done := fastq.ReadFqBig(file); !done; read, done = fastq.ReadFqBig(file) {
		fqs = append(fqs, read)
	}
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}

	//t.ResetTimer()
	//for bench := 0; bench < b.N; bench++ {

	memory := NewMatrixPool()

	var start, stop time.Time

	for read := 0; read < len(fqs); read++ {
		start = time.Now()

		GraphSmithWatermanToGiraf(genome, fqs[read], tiles, tileSize, stepSize, memory, align.HumanChimpTwoScoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)

		stop = time.Now()

		//b.Logf("%s\n", giraf.ToString(result))
	}

	duration := stop.Sub(start)
	// b.ReportAllocs()
	b.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
	//}
	os.Remove(fqOne)

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}

// func TestRightAlignGraphTraversal(t *testing.T) {
// 	// Setup a simple genome graph for testing
// 	genome := EmptyGraph()

// 	var n0, n1, n2 Node
// 	var e0, e1 Edge

// 	// Graph nodes: 	AAAAAT -> ATCG -> CGTTTTTT
// 	// Input Read: 		ATATCGCGT
// 	// Output Expected:	AT->ATCG->CGT
// 	n0 = Node{
// 		Id:  0,
// 		Seq: dna.StringToBases("AAAAACATGTTTT"),
// 	}

// 	n1 = Node{
// 		Id:  1,
// 		Seq: dna.StringToBases("CAT"),
// 	}

// 	n2 = Node{
// 		Id:  2,
// 		Seq: dna.StringToBases("GTTTT"),
// 	}

// 	// Make Edges
// 	e0 = Edge{
// 		Dest: &n1,
// 		Prob: 1,
// 	}

// 	e1 = Edge{
// 		Dest: &n2,
// 		Prob: 1,
// 	}

// 	// Define Paths
// 	n0.Next = append(n0.Next, e0)
// 	n1.Next = append(n1.Next, e1)
// 	n1.Prev = append(n1.Prev, Edge{Dest: &n0, Prob: 1})
// 	n2.Prev = append(n2.Prev, Edge{Dest: &n1, Prob: 1})

// 	genome.Nodes = append(genome.Nodes, n0, n1, n2)

// 	// Define the rest of the test setup...
// 	pool := &sync.Pool{
// 		New: func() interface{} {
// 			return &dnaPool{}
// 		},
// 	}

// 	read := dna.StringToBases("ATGTT")
// 	matrix := NewSwMatrix(defaultMatrixSize)

// 	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 6, Op: 'M'}}
// 	expectedRightTargetStart := 3
// 	expectedRightQueryStart := 6
// 	expectedRightPath := []uint32{1,2}

// 	scoreMatrix := align.HumanChimpTwoScoreMatrix // Assuming this is defined elsewhere

// 	// Perform right traversal starting from the first node (n0)
// 	rightAlign, rightScore, rightTargetStart, rightQueryStart, rightPath := LocalDynamicAln(&n1, read, 0, []uint32{}, len(read), read, scoreMatrix, &matrix, scoreKeeper{},  dynamicScoreKeeper{}, pool)

// 	// Expected results for RightDirection
// 	expectedRightScore := int64(560) // Simplified expected score

// 	// Assert RightDirection results
// 	if rightScore != expectedRightScore {
// 		t.Errorf("RightDirection: expected score %d, got %d", expectedRightScore, rightScore)
// 	}
// 	if !reflect.DeepEqual(cigar.ByteCigarToString(rightAlign), cigar.ByteCigarToString(expectedRightAlignment)) {
// 		t.Errorf("RightDirection: expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
// 	}
// 	if rightTargetStart != expectedRightTargetStart || rightQueryStart != expectedRightQueryStart {
// 		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetStart, expectedRightQueryStart, rightTargetStart, rightQueryStart)
// 	}
// 	if PathToString(rightPath) != PathToString(expectedRightPath) {
// 		t.Errorf("Expected path %+v, got %+v", expectedRightPath, rightPath)
// 	}
// }

// func TestNodeSeqTraversal(t *testing.T) {
// 	// Create a test node
// 	node := Node{
// 		Id:  0,
// 		Seq: dna.StringToBases("ACGT"),
// 	}

// 	// Test leftTraversal
// 	extension := 2
// 	position := 2
// 	seq := dna.StringToBases("TT")
// 	direction := leftTraversal
// 	expectedResult := dna.StringToBases("ACTT")
// 	result := nodeSeqTraversal(&node, extension, position, seq, []dna.Base{}, direction)
// 	if !reflect.DeepEqual(result, expectedResult) {
// 		t.Errorf("LeftTraversal: expected %v, got %v", expectedResult, result)
// 	}
// 	// Test rightTraversal
// 	extension = 2
// 	position = 1
// 	seq = dna.StringToBases("TT")
// 	direction = rightTraversal
// 	expectedResult = dna.StringToBases("TTCG")
// 	result = nodeSeqTraversal(&node, extension, position, seq, []dna.Base{}, direction)
// 	if !reflect.DeepEqual(result, expectedResult) {
// 		t.Errorf("RightTraversal: expected %v, got %v", expectedResult, result)
// 	}
// }
