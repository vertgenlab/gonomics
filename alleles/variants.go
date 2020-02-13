package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
	"runtime"
	"strings"
	"sync"
)

type BatchAlleleCount struct {
	Sample 	string
	Ref    	dna.Base
	Counts 	int32
	BaseA  	[]int32
	BaseC  	[]int32
	BaseG  	[]int32
	BaseT  	[]int32
	Indel  	[]Indel
}

type ScoreInput struct {
	a             int32
	b             int32
	c             int32
	d             int32
	loc           Location
	inStruct      *BatchAlleleCount
	altbase       string
	indelslicepos int
	afThreshold   float64
	sigThreshold  float64
}

// Map structure: map[Chromosome]map[Position][Sample]*BatchAlleleCount
type BatchSampleMap map[Location][]*BatchAlleleCount

// Input a directory filled with AlleleCount (.ac) files. Merges the files into a nested map (chr->pos->[sample]data)
func CreateBatchSampleMap(inDirectory string) BatchSampleMap {

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)
	var SampleName string
	var current *BatchAlleleCount
	var fileCount int = 1
	var i, j int

	SampleMap := make(BatchSampleMap)

	for _, file := range files {

		SampleName = file.Name()
		SamplePath := fmt.Sprintf("%s%s", inDirectory, SampleName)
		log.Printf("#\n# File %d\n", fileCount)
		AlleleCounts := ReadVcfToAlleleCounts(SamplePath)
		fileCount++

		// Loop through input map and deposit info into output map
		for loc, alleles := range AlleleCounts {

			// Check if position is in the map, if not then add the background struct
			_, ok := SampleMap[Location{loc.Chr, loc.Pos}]
			if !ok {
				current = &BatchAlleleCount{"Background", alleles.Ref, 0, make([]int32, 3), make([]int32, 3), make([]int32, 3), make([]int32, 3), make([]Indel, 1)}
				current.Indel[0].Count = make([]int32, 1)
				SampleMap[Location{loc.Chr, loc.Pos}] = append(SampleMap[Location{loc.Chr, loc.Pos}], current)
			}

			//add in and increment background as entry zero in batch allele count slice
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseA[0] = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseA[0] + alleles.BaseA[0]
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseC[0] = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseC[0] + alleles.BaseC[0]
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseG[0] = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseG[0] + alleles.BaseG[0]
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseT[0] = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseT[0] + alleles.BaseT[0]
			SampleMap[Location{loc.Chr, loc.Pos}][0].Counts = SampleMap[Location{loc.Chr, loc.Pos}][0].Counts + alleles.Counts

			// Loop through each indel and see if it is already in the background struct, if not then append it
			for i = 0; i < len(alleles.Indel); i++ {
				var Match bool = false
				for j = 0; j < len(SampleMap[Location{loc.Chr, loc.Pos}][0].Indel); j++ {
					if dna.CompareSeqsIgnoreCase(alleles.Indel[i].Alt, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles.Indel[i].Ref, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Ref) == 0 {
						SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Count[0] = SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Count[0] + alleles.Indel[i].Count[0]
						Match = true
					}
				}
				if Match == false {
					SampleMap[Location{loc.Chr, loc.Pos}][0].Indel = append(SampleMap[Location{loc.Chr, loc.Pos}][0].Indel, alleles.Indel[i])
				}
			}

			current = &BatchAlleleCount{
				Sample: SampleName,
				Ref:    alleles.Ref,
				Counts: alleles.Counts,
				BaseA:  alleles.BaseA,
				BaseC:  alleles.BaseC,
				BaseG:  alleles.BaseG,
				BaseT:  alleles.BaseT,
				Indel:  alleles.Indel}

			SampleMap[Location{loc.Chr, loc.Pos}] = append(SampleMap[Location{loc.Chr, loc.Pos}], current)
		}
	}

	return SampleMap
}

// Call Calculate pValues with fisher's exact test and export in a vcf.Vcf struct with pValue stored in Qual field. Filters on sigThreshold.
func ScoreVariants(input BatchSampleMap, sigThreshold float64, afThreshold float64, numGoRoutines int) []*vcf.Vcf {

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
				loc:           loc,
				sigThreshold:  sigThreshold,
				indelslicepos: 0}

			if alleles[i].Ref != dna.A && passStrandBias(alleles[i].BaseA[1], alleles[i].BaseA[2]) {
				fetInput.c = cA
				fetInput.d = dA
				fetInput.altbase = "A"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.C && passStrandBias(alleles[i].BaseA[1], alleles[i].BaseA[2]) {
				fetInput.c = cC
				fetInput.d = dC
				fetInput.altbase = "C"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.G && passStrandBias(alleles[i].BaseA[1], alleles[i].BaseA[2]) {
				fetInput.c = cG
				fetInput.d = dG
				fetInput.altbase = "G"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.T && passStrandBias(alleles[i].BaseA[1], alleles[i].BaseA[2]) {
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

				if passStrandBias(alleles[i].Indel[j].Count[1], alleles[i].Indel[j].Count[2]) == false {
					continue
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
		fmt.Println(runtime.NumGoroutine())
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

// Calculates strand bias
func passStrandBias(alpha int32, beta int32) bool {
	val := float64(alpha) / float64(alpha+beta)
	if val < 0.95 && val > 0.05 {
		return true
	} else {
		return false
	}
}

// Includes logic to exclude putative variants for which fishers exact test is unnecessary (e.g. alt allele count = 0) and exports as vcf.Vcf

func score(input ScoreInput) *vcf.Vcf{
	var p float64
	var answer *vcf.Vcf

	switch {
	// If alternate allele is zero then there is no variant and score is 1
	case input.c == 0:
		p = 1

	// If a = b and c = d then it is testing itself and should return 1
	case input.a == input.b && input.c == input.d:
		p = 1

	// If the allele frequency of d > c then p is 1
	case float64(input.c)/float64(input.c+input.a) < float64(input.d)/float64(input.d+input.b):
		p = 1

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be excluded
	case float64(input.c)/float64(input.c+input.a) < input.afThreshold:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(int(input.a), int(input.b), int(input.c), int(input.d), true)
	}

	if p < input.sigThreshold {

		var sampleField []string = make([]string, 0)
		sampleField = append(sampleField, fmt.Sprintf("%s:%d:%d:%d", input.inStruct.Sample, input.a, input.c, input.inStruct.Counts))
		answer = &vcf.Vcf{
			Chr:     input.loc.Chr,
			Id:      ".",
			Qual:    p,
			Filter:  ".",
			Format:  "Sample:RefCount:AltCount:Cov",
			Sample: sampleField}

		switch input.altbase {
		case "A", "C", "G", "T":
			answer.Pos = input.loc.Pos + 1
			answer.Info = "."
			answer.Ref = dna.BaseToString(input.inStruct.Ref)
			answer.Alt = input.altbase

		case "Indel":
			answer.Pos = input.loc.Pos
			answer.Info = "."
			answer.Ref = dna.BasesToString(input.inStruct.Indel[input.indelslicepos].Ref)
			answer.Alt = dna.BasesToString(input.inStruct.Indel[input.indelslicepos].Alt)
		}
	}
	return answer
}

// Inputs a vcf file annotated by SnpEff and exports as a csv
func EffToCSV(inFile string, outFile string) {
	data := vcf.Read(inFile)

	// Write header to output file
	output := fileio.MustCreate(outFile)
	defer output.Close()
	fmt.Fprintf(output, "Sample,Chr,Pos,Ref,Alt,pVal,Gene,DNA,Protein,RefCount,AltCount,Coverage,Consequence,Prediction,Transcript\n")

	var i int
	for i = 0; i < len(data); i++ {

		Sample := strings.Split(data[i].Sample[0], ":")
		Info := strings.Split(data[i].Info, "=")

		if len(Info) < 2 {
			continue
		}

		// Eff Example Line: ANN=A|missense_variant|MODERATE|TGFBR3|TGFBR3|transcript|NM_003243.4|protein_coding|13/17|c.2050G>T|p.Asp684Tyr|2565/6467|2050/2556|684/851||
		Ann := strings.Split(Info[1], "|")
		Cov := strings.Split(Sample[3], "\n")

		fmt.Fprintf(output, "%s,%s,%d,%s,%s,%.3v,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
			Sample[0],    // Sample
			data[i].Chr,  // Chr
			data[i].Pos,  // Pos
			data[i].Ref,  // Ref
			data[i].Alt,  // Alt
			data[i].Qual, // pVal
			Ann[3],       // Gene
			Ann[9],       // DNA
			Ann[10],      // Protein
			Sample[1],    // RefCount
			Sample[2],    // AltCount
			Cov[0],       // Coverage
			Ann[1],       // Consequence
			Ann[2],       // Prediction
			Ann[6])       // Transcript
	}
}
