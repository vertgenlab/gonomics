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
		Sample: []string{}}

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

type GraphLocation struct {
	Node 	*simpleGraph.Node
	Pos 	int64
}

func GraphCountAlleles(graph *simpleGraph.SimpleGraph, samFilename string, minMapQ int64) GraphSampleMap {
	AlleleMap := make(GraphSampleMap)

	var i, k int32
	var j int
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64
	var Match bool

	sam.ReadHeader(samFile)

	var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		readPath := StringToPath(aln.Extra)

		if progressMeter%10000 == 0 {
			log.Printf("#%s: Read %d Alignments\n", samFilename, progressMeter)
		}
		progressMeter++

		// If read is unmapped then go to the next alignment
		if aln.Cigar[0].Op == '*' {
			continue
		}

		// If mapping quality is less than the threshold then go to next alignment
		if aln.MapQ < minMapQ {
			continue
		}

		SeqIndex = 0
		RefIndex = aln.Pos - 1

		currNode := 0
		ref := graph.Nodes[readPath[currNode]]

		for i = 0; i < int32(len(aln.Cigar)); i++ {
			currentSeq = aln.Seq

			//Handle deletion relative to ref
			//Each position deleted is annotated with counts + 1
			if aln.Cigar[i].Op == 'D' {
				OrigRefIndex = RefIndex
				OrigNode := ref
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				if OrigRefIndex == 0 {
					prevNodeSeq := graph.Nodes[readPath[currNode - 1]].Seq
					indelSeq[0] = prevNodeSeq[len(prevNodeSeq) - 1]
				} else {
					indelSeq[0] = OrigNode.Seq[OrigRefIndex - 1]
				}


				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					// If the position has already been added to the map, move along
					_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

					// If the position is NOT in the map, initialize
					if !ok {
						AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
							Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					// Keep track of deleted sequence
					indelSeq = append(indelSeq, ref.Seq[RefIndex])

					AlleleMap[GraphLocation{ref, RefIndex}].Counts++

					if RefIndex + 1 == int64(len(ref.Seq)) {
						RefIndex = 0
						currNode++
						if currNode + 1 <= len(readPath) {
							ref = graph.Nodes[readPath[currNode]]
						} else {break}
					} else {
						RefIndex++
					}

				}

				Match = false
				for j = 0; j < len(AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel); j++ {
					// If the deletion has already been seen before, increment the existing entry
					// For a deletion the indelSeq should match the Ref
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Ref) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Alt) == 0 {
						AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[2]++
						}

						Match = true
						break
					}
				}

				// If the deletion has not been seen before, then append it to the Del slice
				// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
				if Match == false {

					currentIndel = Indel{indelSeq, indelSeq[:1], make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == false {
						currentIndel.Count[2]++
					}
					AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel = append(AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel, currentIndel)
				}

				//Handle insertion relative to ref
				//The base after the inserted sequence is annotated with an Ins read
			} else if aln.Cigar[i].Op == 'I' {

				// If the position has already been added to the map, move along
				_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

				// If the position is NOT in the map, initialize
				if !ok {
					AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
						Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
				}

				// Loop through read sequence and keep track of the inserted bases
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				if RefIndex == 0 {
					prevNodeSeq := graph.Nodes[readPath[currNode - 1]].Seq
					indelSeq[0] = prevNodeSeq[len(prevNodeSeq) - 1]
				} else {
					indelSeq[0] = ref.Seq[RefIndex - 1]
				}

				// Keep track of inserted sequence by moving along the read
				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
					indelSeq = append(indelSeq, currentSeq[SeqIndex])
					SeqIndex++
				}

				Match = false
				for j = 0; j < len(AlleleMap[GraphLocation{ref, RefIndex}].Indel); j++ {
					// If the inserted sequence matches a previously inserted sequence, then increment the count
					// For an insertion, the indelSeq should match the Alt
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Ref) == 0 {
						AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[2]++
						}
						Match = true
						break
					}
				}

				if Match == false {
					currentIndel = Indel{indelSeq[:1], indelSeq, make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == true {
						currentIndel.Count[2]++
					}
					AlleleMap[GraphLocation{ref, RefIndex}].Indel = append(AlleleMap[GraphLocation{ref, RefIndex}].Indel, currentIndel)
				}

				// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

				//Handle matching pos relative to ref
			} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					//if the position has already been added to the matrix, move along
					_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

					//if the position is NOT in the matrix, add it
					if !ok {
						AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
							Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					switch currentSeq[SeqIndex] {
					case dna.A:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseA[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseA[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseA[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.T:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseT[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseT[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseT[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.G:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseG[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseG[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseG[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.C:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseC[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseC[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseC[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					}
					SeqIndex++
					if RefIndex + 1 == int64(len(ref.Seq)) {
						RefIndex = 0
						currNode++
						if currNode + 1 <= len(readPath) {
							ref = graph.Nodes[readPath[currNode]]
						} else {break}
					} else {
						RefIndex++
					}
				}
			} else if aln.Cigar[i].Op != 'H' {
				SeqIndex = SeqIndex + aln.Cigar[i].RunLength
			}
		}
	}


	return AlleleMap
}

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
				write.Indel[j].Count = addSlice(write.Indel[j].Count, read.Indel[i].Count)
				Match = true
				break
			}
		}

		if Match == false {
			write.Indel = append(write.Indel, read.Indel[i])
		}
	}
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
				BaseA: 		make([]int32, 3),
				BaseC: 		make([]int32, 3),
				BaseG: 		make([]int32, 3),
				BaseT: 		make([]int32, 3),
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
				a[i-1] = alleles[i].BaseA[0]
				b[i-1] = alleles[0].BaseA[0] - alleles[i].BaseA[0]
			}
		case dna.C:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseC[0]
				b[i-1] = alleles[0].BaseC[0] - alleles[i].BaseC[0]
			}
		case dna.G:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseG[0]
				b[i-1] = alleles[0].BaseG[0] - alleles[i].BaseG[0]
			}
		case dna.T:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseT[0]
				b[i-1] = alleles[0].BaseT[0] - alleles[i].BaseT[0]
			}
		default:
			continue
		}

		// Loop through samples and generate scores
		// For loop starts at index 1 because index zero is the background values
		for i = 1; i < len(alleles); i++ {

			// Retrieve Values for c
			cA = alleles[i].BaseA[0]
			cC = alleles[i].BaseC[0]
			cG = alleles[i].BaseG[0]
			cT = alleles[i].BaseT[0]

			// Retrieve Values for d
			dA = alleles[0].BaseA[0] - alleles[i].BaseA[0]
			dC = alleles[0].BaseC[0] - alleles[i].BaseC[0]
			dG = alleles[0].BaseG[0] - alleles[i].BaseG[0]
			dT = alleles[0].BaseT[0] - alleles[i].BaseT[0]

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
				doesPassStrandBias = passStrandBias(alleles[i].BaseA[1], alleles[i].BaseA[2])
			}

			if alleles[i].Ref != dna.A && doesPassStrandBias {
				fetInput.c = cA
				fetInput.d = dA
				fetInput.altbase = "A"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseC[1], alleles[i].BaseC[2])
			}

			if alleles[i].Ref != dna.C && doesPassStrandBias {
				fetInput.c = cC
				fetInput.d = dC
				fetInput.altbase = "C"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseG[1], alleles[i].BaseG[2])
			}

			if alleles[i].Ref != dna.G && doesPassStrandBias {
				fetInput.c = cG
				fetInput.d = dG
				fetInput.altbase = "G"
				inputChan <- fetInput
			}

			if paired == true {
				doesPassStrandBias = passStrandBias(alleles[i].BaseT[1], alleles[i].BaseT[2])
			}

			if alleles[i].Ref != dna.T && doesPassStrandBias {
				fetInput.c = cT
				fetInput.d = dT
				fetInput.altbase = "T"
				inputChan <- fetInput
			}

			// Calculate p for each Indel
			for j = 0; j < len(alleles[i].Indel); j++ {
				cIndel = alleles[i].Indel[j].Count[0]
				// Find Indel in the background Indel slice
				for l = 0; l < len(alleles[0].Indel); l++ {
					if dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Alt, alleles[0].Indel[l].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Ref, alleles[0].Indel[l].Ref) == 0 {
						dIndel = alleles[0].Indel[l].Count[0] - alleles[i].Indel[j].Count[0]
						break
					}
				}

				if paired == true {
					if passStrandBias(alleles[i].Indel[j].Count[1], alleles[i].Indel[j].Count[2]) == false {
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
