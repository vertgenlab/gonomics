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
	BaseAF  int32
	BaseCF  int32
	BaseGF  int32
	BaseTF  int32
	BaseAR  int32
	BaseCR  int32
	BaseGR  int32
	BaseTR  int32
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
				current = &BatchAlleleCount{"Background", alleles.Ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, make([]Indel, 1)}
				current.Indel[0].CountF = 0
				current.Indel[0].CountR = 0
				SampleMap[Location{loc.Chr, loc.Pos}] = append(SampleMap[Location{loc.Chr, loc.Pos}], current)
			}

			//add in and increment background as entry zero in batch allele count slice
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseAF = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseAF + alleles.BaseAF
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseCF = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseCF + alleles.BaseCF
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseGF = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseGF + alleles.BaseGF
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseTF = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseTF + alleles.BaseTF
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseAR = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseAR + alleles.BaseAR
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseCR = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseCR + alleles.BaseCR
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseGR = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseGR + alleles.BaseGR
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseTR = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseTR + alleles.BaseTR
			SampleMap[Location{loc.Chr, loc.Pos}][0].Counts = SampleMap[Location{loc.Chr, loc.Pos}][0].Counts + alleles.Counts

			// Loop through each indel and see if it is already in the background struct, if not then append it
			for i = 0; i < len(alleles.Indel); i++ {
				var Match bool = false
				for j = 0; j < len(SampleMap[Location{loc.Chr, loc.Pos}][0].Indel); j++ {
					if dna.CompareSeqsIgnoreCase(alleles.Indel[i].Alt, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles.Indel[i].Ref, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Ref) == 0 {
						SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].CountF = SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].CountF + alleles.Indel[i].CountF
						SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].CountR = SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].CountR + alleles.Indel[i].CountR
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
				BaseAF:  alleles.BaseAF,
				BaseCF:  alleles.BaseCF,
				BaseGF:  alleles.BaseGF,
				BaseTF:  alleles.BaseTF,
				BaseAR:  alleles.BaseAR,
				BaseCR:  alleles.BaseCR,
				BaseGR:  alleles.BaseGR,
				BaseTR:  alleles.BaseTR,
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
				a[i-1] = alleles[i].BaseAF + alleles[i].BaseAR
				b[i-1] = (alleles[0].BaseAF + alleles[0].BaseAR) - (alleles[i].BaseAF + alleles[i].BaseAR)
			}
		case dna.C:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseCF + alleles[i].BaseCR
				b[i-1] = (alleles[0].BaseCF + alleles[0].BaseCR) - (alleles[i].BaseCF + alleles[i].BaseCR)
			}
		case dna.G:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseGF + alleles[i].BaseGR
				b[i-1] = (alleles[0].BaseGF + alleles[0].BaseGR) - (alleles[i].BaseGF + alleles[i].BaseGR)
			}
		case dna.T:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseTF + alleles[i].BaseTR
				b[i-1] = (alleles[0].BaseTF + alleles[0].BaseTR) - (alleles[i].BaseTF + alleles[i].BaseTR)
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
			dA = (alleles[0].BaseAF + alleles[0].BaseAR) - (alleles[i].BaseAF + alleles[i].BaseAR)
			dC = (alleles[0].BaseCF + alleles[0].BaseCR) - (alleles[i].BaseCF + alleles[i].BaseCR)
			dG = (alleles[0].BaseGF + alleles[0].BaseGR) - (alleles[i].BaseGF + alleles[i].BaseGR)
			dT = (alleles[0].BaseTF + alleles[0].BaseTR) - (alleles[i].BaseTF + alleles[i].BaseTR)

			// Generate Scores

			fetInput := ScoreInput{
				a:             a[i-1],
				b:             b[i-1],
				afThreshold:   afThreshold,
				inStruct:      alleles[i],
				loc:           loc,
				sigThreshold:  sigThreshold,
				indelslicepos: 0}

			if alleles[i].Ref != dna.A && passStrandBias(alleles[i].BaseAF, alleles[i].BaseAR) {
				fetInput.c = cA
				fetInput.d = dA
				fetInput.altbase = "A"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.C && passStrandBias(alleles[i].BaseCF, alleles[i].BaseCR) {
				fetInput.c = cC
				fetInput.d = dC
				fetInput.altbase = "C"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.G && passStrandBias(alleles[i].BaseGF, alleles[i].BaseGR) {
				fetInput.c = cG
				fetInput.d = dG
				fetInput.altbase = "G"
				inputChan <- fetInput
			}
			if alleles[i].Ref != dna.T && passStrandBias(alleles[i].BaseTF, alleles[i].BaseTR) {
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
						dIndel = (alleles[0].Indel[l].CountF + alleles[0].Indel[l].CountR) - (alleles[i].Indel[j].CountF + alleles[i].Indel[j].CountR)
						break
					}
				}

				if passStrandBias(alleles[i].Indel[j].CountF, alleles[i].Indel[j].CountR) == false {
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
			Notes: sampleField[0]}

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
	vcfchan := vcf.ReadToChan(inFile)

	// Write header to output file
	output := fileio.MustCreate(outFile)
	defer output.Close()
	fmt.Fprintf(output, "Sample\tChr\tPos\tRef\tAlt\tpVal\tGene\tDNA\tProtein\tRefCount\tAltCount\tCoverage\tConsequence\tPrediction\tTranscript\tTLOD Score\tPopulation AF\tCOSMIC ID\tSIFT pred\tPOLYPHEN pred\tPROVEAN pred\n")

	for data := range(vcfchan) {
		Sample := strings.Split(data.Notes, ":")
		Info := strings.Split(data.Info, ";")
		Chr := strings.Split(data.Chr, ":")

		if len(Info) < 2 {
			continue
		}

		var ANNOTATIONS string = ""
		var TLOD string = ""
		var COSMIC string = ""
		var SIFT string = ""
		var POLYPHEN string = ""
		var POPAF string = ""
		var PROVEAN string = ""

		var testInfo []string

		var j int
		for j = 0; j < len(Info); j++ {
			testInfo = strings.Split(Info[j], "=")
			if testInfo[0] == "ANN" {
				ANNOTATIONS = testInfo[1]
			}

			if testInfo[0] == "dbNSFP_1000Gp1_AF" {
				tmp := strings.Split(testInfo[1], ",")
				POPAF = tmp[0]
			}

			if testInfo[0] == "TLOD" {
				tmp := strings.Split(testInfo[1], ",")
				TLOD = tmp[0]
			}

			if testInfo[0] == "dbNSFP_SIFT_pred" {
				tmp := strings.Split(testInfo[1], ",")
				SIFT = tmp[0]
			}

			if testInfo[0] == "dbNSFP_Polyphen2_HVAR_pred" {
				tmp := strings.Split(testInfo[1], ",")
				POLYPHEN = tmp[0]
			}

			if testInfo[0] == "dbNSFP_COSMIC_ID" {
				tmp := strings.Split(testInfo[1], ",")
				COSMIC = tmp[0]
			}

			if testInfo[0] == "dbNSFP_PROVEAN_pred" {
				tmp := strings.Split(testInfo[1], ",")
				PROVEAN = tmp[0]
			}
		}

		// Eff Example Line: ANN=A|missense_variant|MODERATE|TGFBR3|TGFBR3|transcript|NM_003243.4|protein_coding|13/17|c.2050G>T|p.Asp684Tyr|2565/6467|2050/2556|684/851||
		Ann := strings.Split(ANNOTATIONS, "|")
		Cov := strings.Split(Sample[3], "\n")
		counts := strings.Split(Sample[1], ",")

		fmt.Fprintf(output, "%s\t%s\t%d\t%s\t%s\t%.3v\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			Chr[0],    // Sample
			Chr[1],  // Chr
			data.Pos,  // Pos
			data.Ref,  // Ref
			data.Alt,  // Alt
			data.Qual, // pVal
			Ann[3],       // Gene
			Ann[9],       // DNA
			Ann[10],      // Protein
			counts[0],    // RefCount
			counts[1],    // AltCount
			Cov[0],       // Coverage
			Ann[1],       // Consequence
			Ann[2],       // Prediction
			Ann[6],      // Transcript
			TLOD,
			POPAF,
			COSMIC,
			SIFT,
			POLYPHEN,
			PROVEAN)
	}
}


