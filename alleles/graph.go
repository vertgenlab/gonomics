package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
	"strconv"
	"strings"
	"sync"
)

// Functions to find vars that exist in graph
type SampleData struct {
	Sample			string
	Counts 			[]uint32
	GenomeCount 	uint32
}

type BatchData struct {
	Samples 		[]string
	Counts 			[][]uint32
	GenomeCounts 	[]uint32
}

type NodeP struct {
	Sample 	string
	Node 	*simpleGraph.Node
	pVal	float64
}

func StringToPath(input string) []uint32 {
	answer := make([]uint32, 0)
	words := strings.Split(input, ":")
	for i := 0; i < len(words); i++ {
		node, _ := strconv.ParseUint(words[i], 10, 32)
		answer = append(answer, uint32(node))
	}
	return answer
}

func addPathToCount(data []uint32, path []uint32) {
	for i := 0; i < len(path); i++ {
		data[path[i]]++
	}
}

func CountPaths(graph *simpleGraph.SimpleGraph, samFilename string, minMapQ int64) ([]uint32, uint32) {
	answer := make([]uint32, len(graph.Nodes))
	var done = false
	var aln *sam.SamAln

	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	sam.ReadHeader(samFile)

	var GenomeCounts uint32 = 0

	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		if aln.MapQ < minMapQ {continue}
		GenomeCounts++
		addPathToCount(answer, StringToPath(aln.Extra))
	}
	return answer, GenomeCounts
}

func CountPathsInDir(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64) *BatchData {
	var wg sync.WaitGroup
	receivePathCount := make(chan *SampleData)

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)

	for _, file := range files {
		filePath := fmt.Sprintf("%s%s", inDirectory, file.Name())
		// Count Paths for Each Sample in Directory
		wg.Add(1)
		go func(filePath string, fileName string) {
			pathCount, genomeCount := CountPaths(graph, filePath, minMapQ)
			tmp := &SampleData{fileName, pathCount, genomeCount}
			receivePathCount <- tmp
			wg.Done()
		}(filePath, file.Name())
	}

	go func() {
		wg.Wait()
		close(receivePathCount)
	}()

	// Get path counts for each sample
	pathCounts := make([][]uint32, 0)
	samples := make([]string, 0)
	genomeCounts := make([]uint32, 0)
	for count := range receivePathCount {
		samples = append(samples, count.Sample)
		pathCounts = append(pathCounts, count.Counts)
		genomeCounts = append(genomeCounts, count.GenomeCount)
	}

	answer := &BatchData{samples, pathCounts, genomeCounts}
	return answer
}

func calcRareNode(wg *sync.WaitGroup, sendResult chan *NodeP, start *simpleGraph.Node, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) {
	defer wg.Done()
	pathCounts := sampleData.Counts
	genomeCounts := sampleData.GenomeCounts

	// If there were no reads in any adjacent nodes then exit
	// TODO: find some way to determine population frequency and test vs threshold
	//if majorLocalNode == nil {
	//	return
	//}

	// If the prob of entry from the majorLocalNode is > maxPopFreq then exit.
	// No probabilty for first node in graph


	//var prob float64
	//if len(start.Prev) > 0 {
	//	prob = majorPrevEdge.Prob
	//} else {
	//	prob = 1
	//}
	//
	//if prob > maxPopFreq {
	//	return
	//}


	// Begin testing on each sample
	// test is for the matrix:
	// [a b]
	// [c d]
	// a = Samples Genome Count
	// b = Background Genome Count - Samples Genome Count
	// c = Samples Alt Path Count
	// d = Background Alt Path Count - Samples Alt Path Count
	for i := 0; i < len(pathCounts); i++ {
		a := int(genomeCounts[i])
		c := int(pathCounts[i][start.Id])
		var b, d int

		for j := 0; j < len(genomeCounts); j++ {
			b += int(genomeCounts[i])

		}

		for k := 0; k < len(pathCounts); k++ {
			d += int(pathCounts[k][start.Id])
		}

		b = b - a
		d = d - c

		p := numbers.FisherExact(a, b, c, d, true)
		switch {
			// If alternate allele is zero then there is no variant and should exit
			case c == 0:
				return

			// If a = b and c = d then it is testing itself and should exit
			case a == b && c == d:
				return
		}

		if p < minPval {
			answer := &NodeP{sampleData.Samples[i], start, p}
			sendResult <- answer
		}
	}
}

func FindRareNodes(graph *simpleGraph.SimpleGraph, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) []*NodeP {
	var answer = make([]*NodeP, 0)
	result := make(chan *NodeP)
	var wg sync.WaitGroup
	var i, j int

	pathCounts := sampleData.Counts

	// Note: pathCounts is a [sample][nodes]
	// For each node ---->
	for i = 0; i < len(pathCounts[0]); i++ {
		var areThereReads = false
		// For each sample ---->
		for j = 0; j < len(pathCounts); j++ {
			if pathCounts[i][j] != 0 {
				areThereReads = true
				break
			}
		}

		// Ignore node if no reads are mapped to it
		if areThereReads == false {
			continue
		}
		wg.Add(1)
		go calcRareNode(&wg, result, graph.Nodes[i], sampleData, maxPopFreq, minReadFreq, minPval)
	}

	go func() {
		wg.Wait()
		close(result)
	}()

	for node := range result {
		answer = append(answer, node)
	}
	return answer
}

func nodeToVcf(node *NodeP) *vcf.Vcf {
	answer := &vcf.Vcf{
		Chr: 	node.Node.Name,
		Pos: 	0,
		Id: 	"",
		Ref: 	"",
		Alt: 	dna.BasesToString(node.Node.Seq),
		Qual: 	0,
		Filter: "",
		Info: 	fmt.Sprintf("Sample:%s, p:%e", node.Sample, node.pVal),
		Format: "",
		Notes: ""}

	return answer
}

func NodesToVcf(nodes []*NodeP) []*vcf.Vcf {
	answer := make([]*vcf.Vcf, 0)
	var wg sync.WaitGroup
	var receiveVcf = make(chan *vcf.Vcf)

	for i := 0; i < len(nodes); i++ {
		wg.Add(1)
		currNode := nodes[i]
		go func() {
			receiveVcf <- nodeToVcf(currNode)
			wg.Done()
		}()
	}

	go func() {
		wg.Wait()
		close(receiveVcf)
	}()

	for record := range receiveVcf {
		answer = append(answer, record)
	}

	return answer
}

func FindVarsInsideGraph(graph *simpleGraph.SimpleGraph, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) []*vcf.Vcf {
	var answer []*vcf.Vcf

	// Find rare nodes in graph
	rareNodes := FindRareNodes(graph, sampleData, maxPopFreq, minReadFreq, minPval)

	answer = NodesToVcf(rareNodes)

	return answer
}


// Functions to find vars that DO NOT exist in graph
type GraphSampleMap map[GraphLocation]*AlleleCount

type BatchGraphSampleMap map[GraphLocation][]*BatchAlleleCount

func GraphCountAllelesInDir(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64) BatchGraphSampleMap {
	answer := make(BatchGraphSampleMap)

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)

	for _, file := range files {
		filePath := fmt.Sprintf("%s%s", inDirectory, file.Name())
		alleleCount := GraphCountAlleles(graph, filePath, minMapQ)
		appendBGSM(answer, alleleCount, file.Name())
	}

	return answer
}

func batchAddIndels(write *BatchAlleleCount, read *AlleleCount) {
	var Match bool
	var i, j int

	Match = false
	for i = 0; i < len(read.Indel); i++ {
		for j = 0; j < len(write.Indel); j++ {
			// If the deletion has already been seen before, increment the existing entry
			// For a deletion the indelSeq should match the Ref
			if dna.CompareSeqsIgnoreCase(read.Indel[i].Ref, write.Indel[j].Ref) == 0 &&
				dna.CompareSeqsIgnoreCase(read.Indel[i].Alt, write.Indel[j].Alt) == 0 {
				write.Indel[j].CountF = write.Indel[j].CountF + read.Indel[i].CountF
				write.Indel[j].CountR = write.Indel[j].CountR + read.Indel[i].CountR
				Match = true
				break
			}
		}

		if Match == false {
			write.Indel = append(write.Indel, read.Indel[i])
		}
	}
}

func addSlice(a []int32, b []int32) []int32 {
	c := make([]int32, len(a))
	for i := 0; i < len(a); i++ {
		c[i] = a[i] + b[i]
	}
	return c
}

func appendBGSM(OutMap BatchGraphSampleMap, InMap GraphSampleMap, SampleName string) BatchGraphSampleMap {
	for key, value := range InMap {

		_, ok := OutMap[key]
		if !ok {
			OutMap[key] = make([]*BatchAlleleCount, 1)

			bkgd := &BatchAlleleCount{
				Sample:	 	"Background",
				Ref: 		value.Ref,
				Counts: 	0,
				BaseAF: 		0,
				BaseCF: 		0,
				BaseGF: 		0,
				BaseTF: 		0,
				BaseAR: 		0,
				BaseCR: 		0,
				BaseGR: 		0,
				BaseTR: 		0,
				Indel: 		make([]Indel, 1)}


			OutMap[key][0] = bkgd
		}

		OutMap[key][0].Counts += value.Counts
		OutMap[key][0].BaseA = addSlice(OutMap[key][0].BaseA, value.BaseA)
		OutMap[key][0].BaseC = addSlice(OutMap[key][0].BaseC, value.BaseC)
		OutMap[key][0].BaseG = addSlice(OutMap[key][0].BaseG, value.BaseG)
		OutMap[key][0].BaseT = addSlice(OutMap[key][0].BaseT, value.BaseT)
		batchAddIndels(OutMap[key][0], value)

		curr := &BatchAlleleCount{
			Sample:	 	SampleName,
			Ref: 		value.Ref,
			Counts: 	value.Counts,
			BaseA: 		value.BaseA,
			BaseC: 		value.BaseC,
			BaseG: 		value.BaseG,
			BaseT: 		value.BaseT,
			Indel: 		value.Indel}

		OutMap[key] = append(OutMap[key], curr)
	}
	return OutMap
}

func GraphScoreVariants(input BatchGraphSampleMap, sigThreshold float64, afThreshold float64, numGoRoutines int, paired bool) []*vcf.Vcf {

	fmt.Printf("#\n# Calling Variants\n")
	var VariantScores []*vcf.Vcf
	var progressMeter int
	var cA, cC, cG, cT, cIndel, dA, dC, dG, dT, dIndel int32
	var a, b []int32

	// Initialize a buffered channel to send completed vcf structs through
	vcfChannel := make(chan *vcf.Vcf, len(input))
	var wg sync.WaitGroup

	// Start Goroutines
	var threads int
	if len(input) < numGoRoutines {
		threads = len(input)
	} else {
		threads = numGoRoutines
	}

	wg.Add(threads)
	inputChan := make(chan ScoreInput)
	for k := 0; k < threads; k++ {
		go func() {
			for {
				data, ok := <-inputChan
				if !ok {
					wg.Done()
					return
				}
				answer := score(data)
				if answer != nil {
					vcfChannel <- answer
				}
			}
		}()
	}

	// Loop through BatchSampleMap
	for loc, alleles := range input {

		if progressMeter%1000 == 0 {
			log.Printf("# Processed %d Positions\n", progressMeter)
		}
		progressMeter++

		// Must be at least 2 samples in the slice to generate a batch p value.
		// A single sample would be len = 2 because the background values are always element 0
		if len(alleles) <= 2 {
			continue
		}

		// Begin gathering parameters for Fishers Exact Test done in the numbers package
		// test is for the matrix:
		// [a b]
		// [c d]
		// a = Samples Ref Allele Count
		// b = Background Ref Allele Count - Samples Ref Allele Count
		// c = Samples Alt Allele Count
		// d = Background Alt Allele Count - Samples Alt Allele Count

		// Determine Reference Base
		// If ref base is unknown then skip
		var i, j, l int
		a = make([]int32, len(alleles)-1)
		b = make([]int32, len(alleles)-1)
		switch alleles[0].Ref {
		// Loop through samples and gather inputs for a and b
		// For loop starts at index 1 because index zero is the background values
		case dna.A:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseAF + alleles[i].BaseAR
				b[i-1] = alleles[0].BaseAF + alleles[0].BaseAR - a[i-1]
			}
		case dna.C:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseCF + alleles[i].BaseCR
				b[i-1] = alleles[0].BaseCF + alleles[0].BaseCR - a[i-1]
			}
		case dna.G:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseGF + alleles[i].BaseGR
				b[i-1] = alleles[0].BaseGF + alleles[0].BaseGR - a[i-1]
			}
		case dna.T:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseTF + alleles[i].BaseTR
				b[i-1] = alleles[0].BaseTF + alleles[0].BaseTR - a[i-1]
			}
		default:
			continue
		}

		// Loop through samples and generate scores
		// For loop starts at index 1 because index zero is the background values
		for i = 1; i < len(alleles); i++ {

			// Retrieve Values for c
			cA = alleles[i].BaseAF + alleles[i].BaseAR
			cC = alleles[i].BaseCF + alleles[i].BaseCR
			cG = alleles[i].BaseGF + alleles[i].BaseGR
			cT = alleles[i].BaseTF + alleles[i].BaseTR

			// Retrieve Values for d
			dA = alleles[0].BaseAF + alleles[i].BaseAR - cA
			dC = alleles[0].BaseCF + alleles[i].BaseCR - cC
			dG = alleles[0].BaseGF + alleles[i].BaseGR - cG
			dT = alleles[0].BaseTF + alleles[i].BaseTR - cT

			// Generate Scores

			fetInput := ScoreInput{
				a:             a[i-1],
				b:             b[i-1],
				afThreshold:   afThreshold,
				inStruct:      alleles[i],
				loc:           Location{loc.Node.Name, loc.Pos},
				sigThreshold:  sigThreshold,
				indelslicepos: 0}

			var doesPassStrandBias = true

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseAF, alleles[i].BaseAR)
			}

			if alleles[i].Ref != dna.A && doesPassStrandBias {
				fetInput.c = cA
				fetInput.d = dA
				fetInput.altbase = "A"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseCF, alleles[i].BaseCR)
			}

			if alleles[i].Ref != dna.C && doesPassStrandBias {
				fetInput.c = cC
				fetInput.d = dC
				fetInput.altbase = "C"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseGF, alleles[i].BaseGR)
			}

			if alleles[i].Ref != dna.G && doesPassStrandBias {
				fetInput.c = cG
				fetInput.d = dG
				fetInput.altbase = "G"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseTF, alleles[i].BaseTR)
			}

			if alleles[i].Ref != dna.T && doesPassStrandBias {
				fetInput.c = cT
				fetInput.d = dT
				fetInput.altbase = "T"
				inputChan <- fetInput
			}

			// Calculate p for each Indel
			for j = 0; j < len(alleles[i].Indel); j++ {
				cIndel = alleles[i].Indel[j].CountF + alleles[i].Indel[j].CountR
				// Find Indel in the background Indel slice
				for l = 0; l < len(alleles[0].Indel); l++ {
					if dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Alt, alleles[0].Indel[l].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Ref, alleles[0].Indel[l].Ref) == 0 {
						dIndel = alleles[0].Indel[l].CountF + alleles[i].Indel[j].CountR - cIndel
						break
					}
				}

				if paired == true {
					if passStrandBias(alleles[i].Indel[j].CountF, alleles[i].Indel[j].CountR) == false {
						continue
					}
				}

				fetInput.c = cIndel
				fetInput.d = dIndel
				fetInput.altbase = "Indel"
				fetInput.indelslicepos = j
				inputChan <- fetInput
			}
		}
	}

	fmt.Println("# Waiting for Goroutines to finish")

	// Start GoRoutine to monitor for the finish of the wait group then close the output channel
	go func() {
		wg.Wait()
		fmt.Println("# Goroutines finished")
		close(vcfChannel)
	}()

	close(inputChan)
	wg.Wait()

	for answer := range vcfChannel {
		VariantScores = append(VariantScores, answer)
	}

	return VariantScores
}

func FindVarsOutsideGraph(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64, minPval float64, afThreshold float64, numGoRoutines int, paired bool) []*vcf.Vcf {
	var answer []*vcf.Vcf

	alleleCount := GraphCountAllelesInDir(graph, inDirectory, minMapQ)

	answer = GraphScoreVariants(alleleCount, minPval, afThreshold, numGoRoutines, paired)

	return answer
}


// Complete wrapper function for vars inside and outside graph
func GraphVariants(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64, maxPopFreq float64, minReadFreq float64, minPval float64, afThreshold float64, numGoRoutines int, paired bool) []*vcf.Vcf {

	pathCounts := CountPathsInDir(graph, inDirectory, minMapQ)

	VarsInGraph := FindVarsInsideGraph(graph, pathCounts, maxPopFreq, minReadFreq, minPval)
	VarsOutGraph := FindVarsOutsideGraph(graph, inDirectory, minMapQ, minPval, afThreshold, numGoRoutines, paired)

	answer := append(VarsInGraph, VarsOutGraph...)

	return answer
}
