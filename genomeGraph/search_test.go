package genomeGraph

import (
	"flag"
	"log"
	"os"
	"reflect"
	"runtime"
	"runtime/pprof"
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

	memory := MatrixPoolMemory()
	cache := memory.Get().(*MatrixScore)
	defer memory.Put(cache)

	// Call DynamicAln function
	settings := &GraphSettings{
		scores:     align.DefaultScoreMatrix,
		gapPenalty: -100,
		tileSize:   32,
		stepSize:   32,
		extension:  150,
	}
	var expectedScore int64 = 273
	var expectedTEnd, expectedQEnd int = 9, 6
	expectCigar := []cigar.ByteCigar{{RunLen: 6, Op: 'M'}, {RunLen: 3, Op: 'D'}}

	score, cig, endAlpha, endBeta := RightDynamicAln(alpha, beta, settings, cache)

	t.Logf("Alignment Score: %d\nAlignment: %s\n", score, cigar.ByteCigarToString(cig))

	// Check results
	if score != expectedScore {
		t.Errorf("Error: Score %v != %v\n", score, expectedScore)
	}
	if cigar.ByteCigarToString(cig) != cigar.ByteCigarToString(expectCigar) {
		t.Errorf("Error: CIGAR %+v != %+v\n", cig, expectCigar)
	}
	if endAlpha != expectedTEnd || endBeta != expectedQEnd {
		t.Errorf("Error: Ending indices (%d, %d) != (%d, %d) do not match expected values\n", expectedTEnd, expectedQEnd, endAlpha, endBeta)
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

	memory := MatrixPoolMemory()

	cache := memory.Get().(*MatrixScore)
	defer memory.Put(cache)

	read := dna.StringToBases("ATATCGCGT")

	scoreKeeper := scoreKeeper{}

	// expectedRightPath :=  []uint32{1,2}
	// Expected results for RightDirection
	expectedRightScore := int64(850) // Simplified expected score
	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 9, Op: 'M'}}
	expectedRightTargetEnd := 9
	expectedRightQueryEnd := 9

	settings := &GraphSettings{
		scores:     align.HumanChimpTwoScoreMatrix,
		gapPenalty: -100,
		tileSize:   32,
		extension:  len(read),
	}

	rightAlign, rightScore, rightTargetEnd, rightQueryEnd, _ := RightAlignTraversal(&genome.Nodes[0], read, 0, []uint32{}, read, settings, &scoreKeeper, memory)

	// leftAlign, leftScore, leftTargetStart, leftQueryStart, _ := LeftAlignTraversal(&genome.Nodes[0], read, 0, []uint32{}, read, settings, scoreKeeper, memory)

	if rightScore != expectedRightScore {
		t.Errorf("RightDirection: expected score %d, got %d", expectedRightScore, rightScore)
	}
	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(expectedRightAlignment) {
		t.Errorf("RightDirection: expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
	}
	if expectedRightTargetEnd != rightTargetEnd || rightQueryEnd != expectedRightQueryEnd {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetEnd, rightQueryEnd, rightTargetEnd, rightQueryEnd)
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

	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1
	var readLength int = 150
	var mutations int = 12

	settings := &GraphSettings{
		scores:     align.DefaultScoreMatrix,
		gapPenalty: -600,
		tileSize:   tileSize,
		stepSize:   stepSize,
		extension:  readLength,
	}

	genome := Read("testdata/bigGenome.sg")
	fqOne, fqTwo := "testdata/simReads_R1.fq", "testdata/simReads_R2.fq"

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)

	fastq.WritePair(fqOne, fqTwo, simReads)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	fqs := []fastq.FastqBig{}
	file := fileio.NewByteReader(fqOne)
	for read, done := fastq.ReadFqBig(file); !done; read, done = fastq.ReadFqBig(file) {
		fqs = append(fqs, read)
	}

	memory := MatrixPoolMemory()
	seedPool := NewMemSeedPool()

	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}

	start := time.Now()

	for read := 0; read < len(fqs); read++ {
		GraphSmithWatermanToGiraf(genome, fqs[read], tiles, settings, &seedPool, scorekeeper, seedBuildHelper, memory)
	}

	stop := time.Now()

	duration := stop.Sub(start)
	b.ReportAllocs()

	b.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())

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

func TestRightAlignGraphTraversal(t *testing.T) {
	// Setup a simple genome graph for testing
	genome := EmptyGraph()

	var n0, n1, n2 Node
	var e0, e1 Edge

	// Graph nodes: 	AAAAAT -> ATCG -> CGTTTTTT
	// Input Read: 		ATATCGCGT
	// Output Expected:	AT->ATCG->CGT
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("AAAAACATGTTTT"),
	}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("CAT"),
	}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("GTTTT"),
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
	n0.Next = append(n0.Next, e0)
	n1.Next = append(n1.Next, e1)
	n1.Prev = append(n1.Prev, Edge{Dest: &n0, Prob: 1})
	n2.Prev = append(n2.Prev, Edge{Dest: &n1, Prob: 1})

	genome.Nodes = append(genome.Nodes, n0, n1, n2)

	read := dna.StringToBases("ATGTT")

	memory := MatrixPoolMemory()
	settings := &GraphSettings{
		scores:     align.HumanChimpTwoScoreMatrix,
		gapPenalty: -100,
		tileSize:   32,
		stepSize:   32,
		extension:  150,
	}

	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 5, Op: 'M'}}
	expectedRightTargetStart := 5
	expectedRightQueryStart := 5
	expectedRightPath := []uint32{}

	// Perform right traversal starting from the first node (n0)
	rightAlign, rightScore, rightTargetStart, rightQueryStart, rightPath := RightAlignTraversal(&n1, read, 0, []uint32{}, read, settings, &scoreKeeper{}, memory)

	// Expected results for RightDirection
	expectedRightScore := int64(460) // Simplified expected score

	// Assert RightDirection results
	if rightScore != expectedRightScore {
		t.Errorf("RightDirection: expected score %d, got %d", expectedRightScore, rightScore)
	}
	if !reflect.DeepEqual(cigar.ByteCigarToString(rightAlign), cigar.ByteCigarToString(expectedRightAlignment)) {
		t.Errorf("RightDirection: expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
	}
	if rightTargetStart != expectedRightTargetStart || rightQueryStart != expectedRightQueryStart {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetStart, expectedRightQueryStart, rightTargetStart, rightQueryStart)
	}
	if PathToString(rightPath) != PathToString(expectedRightPath) {
		t.Errorf("Expected path %+v, got %+v", expectedRightPath, rightPath)
	}
}

func TestGetTargetBases(t *testing.T) {
	// Test case 1: direction is left
	n := &Node{
		Seq: dna.StringToBases("ACGT"),
	}
	extension := len(n.Seq)
	position := len(n.Seq) - 1
	direction := left

	expected := "ACG"
	result := dna.BasesToString(getTargetBases(n, extension, position, []dna.Base{}, []dna.Base{}, direction))

	if result != expected {
		t.Errorf("Test case 1 failed. Expected: %v, got: %v", expected, result)
	}

	// Test case 2: direction is right
	n = &Node{
		Seq: dna.StringToBases("ACGT"),
	}
	extension = 5
	position = 0

	direction = right

	expected = "ACGT"
	result = dna.BasesToString(getTargetBases(n, extension, position, []dna.Base{}, []dna.Base{}, direction))

	if result != expected {
		t.Errorf("Test case 2 failed. Expected: %v, got: %v", expected, result)
	}
}
