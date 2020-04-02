package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
	"os"
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
//type GraphSampleMap map[GraphLocation]*AlleleCount

//type BatchGraphSampleMap map[GraphLocation][]*BatchAlleleCount

type BatchAllele struct {
	Sample 	string
	Allele 	*Allele
}

type ProgressLock bool

func GraphCountAllelesInDir(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64) chan []*BatchAllele {

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)

	var currLocation *Location
	var receiveAlleles chan *BatchAllele = make(chan *BatchAllele)
	var sendAlleles chan []*BatchAllele = make(chan []*BatchAllele)

	// Get sam header
	filePath := fmt.Sprintf("%s%s", inDirectory, files[0].Name())
	samFile := fileio.EasyOpen(filePath)
	samHeader := sam.ReadHeader(samFile)

	if len(samHeader.Chroms) == 0 {
		log.Fatalln("ERROR: No sam header detected")
	}

	// Define start position
	currLocation = &Location{samHeader.Chroms[0].Name, 0}

	var wg sync.WaitGroup
	locks := make([]*ProgressLock, 0)

	for _, file := range files {
		var lock ProgressLock = false
		go sendOffAlleles(inDirectory, file, graph, minMapQ, currLocation, receiveAlleles, &wg, &lock)
		locks = append(locks, &lock)
	}

	go listenForAlleles(samHeader, receiveAlleles, sendAlleles, currLocation, &wg, locks)

	return sendAlleles
}

func listenForAlleles(samHeader *sam.SamHeader, receiveAlleles chan *BatchAllele, sendAlleles chan []*BatchAllele, currLocation *Location, wg *sync.WaitGroup, locks []*ProgressLock) {
	var i, k int
	var j int64
	answer := make([]*BatchAllele, 0)

	go func() {
		for allele := range receiveAlleles {
			answer = append(answer, allele)
		}
	}()

	for i = 0; i < len(samHeader.Chroms); i++ {
		currLocation.Chr = samHeader.Chroms[i].Name
		for j = 0; j < samHeader.Chroms[i].Size; j++ {
			currLocation.Pos = j
			wg.Add(1)

			for {
				var readyToProgress bool = true
				// Check to see if all workers are ready for next currLocation
				for k = 0; k < len(locks); k++ {
					if *locks[k] == false {
						readyToProgress = false
						break
					}
				}
				// If all workers are ready, then reset the locks, send the allele, and progress
				if readyToProgress == true {
					for k = 0; k < len(locks); k++ {
						*locks[k] = false
					}

					sendAlleles <- answer
					answer = nil
					wg.Done()
					break
				}
			}
		}
	}
	close(sendAlleles)
	close(receiveAlleles)
}

func sendOffAlleles(inDirectory string, file os.FileInfo, graph *simpleGraph.SimpleGraph, minMapQ int64, currLocation *Location, receiveAlleles chan *BatchAllele, wg *sync.WaitGroup, lock *ProgressLock) {
	filePath := fmt.Sprintf("%s%s", inDirectory, file.Name())
	alleleStream := SamToAlleles(filePath, graph, minMapQ)

	for allele := range alleleStream {
		// Loop forever until the allele is sent
		for {
			// If the allele ready to send is at the currLocaiton, then send it
			if allele.Location.Pos == currLocation.Pos && allele.Location.Chr == currLocation.Chr {
				receiveAlleles <- &BatchAllele{file.Name(), allele}
				break
			} else {
				*lock = true
				wg.Wait()
			}
		}
	}
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
/*
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
*/

func FindMatchingIndel(queryIndel *Indel, subjectSlice []Indel) *Indel {
	var i int
	for i = 0; i < len(subjectSlice); i++ {
		if dna.CompareSeqsIgnoreCase(queryIndel.Alt, subjectSlice[i].Alt) == 0 &&
			dna.CompareSeqsIgnoreCase(queryIndel.Ref, subjectSlice[i].Ref) == 0 {
			return &subjectSlice[i]
		}
	}
	return nil
}

func sumBatchAllele(input []*BatchAllele) *AlleleCount {
	var answer = &AlleleCount{input[0].Allele.Count.Ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, make([]Indel, 0)}
	var count *AlleleCount
	var i, j int
	for i = 0; i < len(input); i++ {
		count = input[i].Allele.Count
		answer.Counts += count.Counts
		answer.BaseAF += count.BaseAF
		answer.BaseAR += count.BaseAR
		answer.BaseCF += count.BaseCF
		answer.BaseCR += count.BaseCR
		answer.BaseGF += count.BaseGF
		answer.BaseGR += count.BaseGR
		answer.BaseTF += count.BaseTF
		answer.BaseTR += count.BaseTR

		for j = 0; j < len(count.Indel); j++ {
			indel := FindMatchingIndel(&count.Indel[j], answer.Indel)
			if indel != nil {
				indel.CountF += count.Indel[j].CountF
				indel.CountR += count.Indel[j].CountR
			} else {
				answer.Indel = append(answer.Indel, count.Indel[j])
			}
		}
	}
	return answer
}

func appendAlleleToVcf(vcf *vcf.Vcf, refBase []dna.Base, altBase []dna.Base, data *BatchAllele, p float64) *vcf.Vcf {

	if vcf.Alt == "" {
		vcf.Alt = dna.BasesToString(altBase)
	} else {
		vcf.Alt = vcf.Alt + "," + dna.BasesToString(altBase)
	}

	// Check if it is a SNP or Indel
	// Format order is "Sample:RefCount:AltCount:Cov:pValue"
	if len(refBase) == 1 && len(altBase) == 1 {
		// Case SNP
		// TODO: fill in fields
	} else {
		// case INDEL
		if len(refBase) != 1 {
			vcf.Ref = vcf.Ref + "," + dna.BasesToString(refBase)
		}
		// TODO: fill in fields
	}
	return vcf
}

func GraphScoreVariant(answer chan *vcf.Vcf, input []*BatchAllele, sigThreshold float64, afThreshold float64, paired bool) {

	fmt.Printf("#\n# Calling Variants\n")
	var cA, cC, cG, cT, cIndel, dA, dC, dG, dT, dIndel int32
	var a, b []int32
	var sumCount *AlleleCount

	// Must be at least 2 samples in the slice to generate a batch p value.
	if len(input) < 2 {
		return
	}

	sumCount = sumBatchAllele(input)

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
	a = make([]int32, len(input))
	b = make([]int32, len(input))
	switch input[0].Allele.Count.Ref {
	// Loop through samples and gather inputs for a and b
	// For loop starts at index 1 because index zero is the background values
	case dna.A:
		for i = 0; i < len(input); i++ {
			a[i] = input[i].Allele.Count.BaseAF + input[i].Allele.Count.BaseAR
			b[i] = sumCount.BaseAF + sumCount.BaseAR - a[i]
		}
	case dna.C:
		for i = 0; i < len(input); i++ {
			a[i] = input[i].Allele.Count.BaseCF + input[i].Allele.Count.BaseCR
			b[i] = sumCount.BaseCF + sumCount.BaseCR - a[i]
		}
	case dna.G:
		for i = 0; i < len(input); i++ {
			a[i] = input[i].Allele.Count.BaseGF + input[i].Allele.Count.BaseGR
			b[i] = sumCount.BaseGF + sumCount.BaseGR - a[i]
		}
	case dna.T:
		for i = 0; i < len(input); i++ {
			a[i] = input[i].Allele.Count.BaseTF + input[i].Allele.Count.BaseTR
			b[i] = sumCount.BaseTF + sumCount.BaseTR - a[i]
		}
	default:
		return
	}

	// Loop through samples and generate scores
	// For loop starts at index 1 because index zero is the background values
	for i = 0; i < len(input); i++ {
		var vcfRecord *vcf.Vcf
		vcfRecord.Chr = input[i].Allele.Location.Chr
		vcfRecord.Pos = input[i].Allele.Location.Pos
		vcfRecord.Ref = dna.BaseToString(input[i].Allele.Count.Ref)
		vcfRecord.Format = "Sample:RefCount:AltCount:Cov:pValue"

		count := input[i].Allele.Count

		// Retrieve Values for c
		cA = count.BaseAF + count.BaseAR
		cC = count.BaseCF + count.BaseCR
		cG = count.BaseGF + count.BaseGR
		cT = count.BaseTF + count.BaseTR

		// Retrieve Values for d
		dA = sumCount.BaseAF + sumCount.BaseAR - cA
		dC = sumCount.BaseCF + sumCount.BaseCR - cC
		dG = sumCount.BaseGF + sumCount.BaseGR - cG
		dT = sumCount.BaseTF + sumCount.BaseTR - cT

		// Generate Scores
		var doesPassStrandBias = true

		if paired == true {
			doesPassStrandBias = passStrandBias(count.BaseAF, count.BaseAR)
		}

		if count.Ref != dna.A && doesPassStrandBias {
			p := ScoreVariant(a[i], b[i], cA, dA, afThreshold)
			if p < sigThreshold {
				appendAlleleToVcf(vcfRecord, []dna.Base{count.Ref}, []dna.Base{dna.A}, input[i], p)
			}
		}

		if paired == true {
			doesPassStrandBias = passStrandBias(count.BaseCF, count.BaseCR)
		}

		if count.Ref != dna.C && doesPassStrandBias {
			p := ScoreVariant(a[i], b[i], cC, dC, afThreshold)
			if p < sigThreshold{
				appendAlleleToVcf(vcfRecord, []dna.Base{count.Ref}, []dna.Base{dna.C}, input[i], p)
			}
		}

		if paired == true {
			doesPassStrandBias = passStrandBias(count.BaseGF, count.BaseGR)
		}

		if count.Ref != dna.G && doesPassStrandBias {
			p := ScoreVariant(a[i], b[i], cG, dG, afThreshold)
			if p < sigThreshold{
				appendAlleleToVcf(vcfRecord, []dna.Base{count.Ref}, []dna.Base{dna.G}, input[i], p)
			}
		}

		if paired == true {
			doesPassStrandBias = passStrandBias(count.BaseTF, count.BaseTR)
		}

		if count.Ref != dna.T && doesPassStrandBias {
			p := ScoreVariant(a[i], b[i], cT, dT, afThreshold)
			if p < sigThreshold{
				appendAlleleToVcf(vcfRecord, []dna.Base{count.Ref}, []dna.Base{dna.T}, input[i], p)
			}
		}

		// Calculate p for each Indel
		for j = 0; j < len(count.Indel); j++ {
			cIndel = count.Indel[j].CountF + count.Indel[j].CountR
			// Find Indel in the background Indel slice
			for l = 0; l < len(sumCount.Indel); l++ {
				if dna.CompareSeqsIgnoreCase(count.Indel[j].Alt, sumCount.Indel[l].Alt) == 0 &&
					dna.CompareSeqsIgnoreCase(count.Indel[j].Ref, sumCount.Indel[l].Ref) == 0 {
					dIndel = sumCount.Indel[l].CountF + count.Indel[j].CountR - cIndel
					break
				}
			}

			if paired == true {
				if passStrandBias(count.Indel[j].CountF, count.Indel[j].CountR) == false {
					continue
				}
			}

			p := ScoreVariant(a[i], b[i], cIndel, dIndel, afThreshold)
			if p < sigThreshold{
				appendAlleleToVcf(vcfRecord, count.Indel[j].Ref, count.Indel[j].Alt, input[i], p)
			}
		}
		answer <- vcfRecord
	}
}

func FindVarsOutsideGraph(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64, minPval float64, afThreshold float64, paired bool) chan *vcf.Vcf {
	answer := make(chan *vcf.Vcf)
	alleleStream := GraphCountAllelesInDir(graph, inDirectory, minMapQ)

	go func() {
		for allele := range alleleStream {
			GraphScoreVariant(answer, allele, minPval, afThreshold, paired)
		}
		close(answer)
	}()

	return answer
}

func ScoreVariant(a int32, b int32, c int32, d int32, afThreshold float64) float64 {
	var p float64

	switch {
	// If alternate allele is zero then there is no variant and score is 1
	case c == 0:
		p = 1

	// If a = b and c = d then it is testing itself and should return 1
	case a == b && c == d:
		p = 1

	// If the allele frequency of d > c then p is 1
	case float64(c)/float64(c+a) < float64(d)/float64(d+b):
		p = 1

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be excluded
	case float64(c)/float64(c+a) < afThreshold:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(int(a), int(b), int(c), int(d), true)
	}
	return p
}


// Complete wrapper function for vars inside and outside graph
func GraphVariants(graph *simpleGraph.SimpleGraph, inDirectory string, minMapQ int64, maxPopFreq float64, minReadFreq float64, minPval float64, afThreshold float64, numGoRoutines int, paired bool) []*vcf.Vcf {
	var answer []*vcf.Vcf
	pathCounts := CountPathsInDir(graph, inDirectory, minMapQ)

	VarsInGraph := FindVarsInsideGraph(graph, pathCounts, maxPopFreq, minReadFreq, minPval)
	VarsOutGraph := FindVarsOutsideGraph(graph, inDirectory, minMapQ, minPval, afThreshold, paired)

	for val := range VarsOutGraph {
		answer = append(VarsInGraph, val)
	}

	return answer
}


