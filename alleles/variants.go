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
// TODO: look into limiting goroutines
type BatchAlleleCount struct {
	Sample 	string
	Ref    	dna.Base
	Counts 	int32
	BaseA  	int32
	BaseC  	int32
	BaseG  	int32
	BaseT  	int32
	Indel  	[]Indel
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
				current = &BatchAlleleCount{"Background", alleles.Ref, 0, 0, 0, 0, 0, make([]Indel, 0)}
				SampleMap[Location{loc.Chr, loc.Pos}] = append(SampleMap[Location{loc.Chr, loc.Pos}], current)
			}

			//add in and increment background as entry zero in batch allele count slice
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseA = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseA + alleles.BaseA
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseC = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseC + alleles.BaseC
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseG = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseG + alleles.BaseG
			SampleMap[Location{loc.Chr, loc.Pos}][0].BaseT = SampleMap[Location{loc.Chr, loc.Pos}][0].BaseT + alleles.BaseT
			SampleMap[Location{loc.Chr, loc.Pos}][0].Counts =SampleMap[Location{loc.Chr, loc.Pos}][0].Counts + alleles.Counts

			// Loop through each indel and see if it is already in the background struct, if not then append it
			for i = 0; i < len(alleles.Indel); i++ {
				var Match bool = false
				for j = 0; j < len(SampleMap[Location{loc.Chr, loc.Pos}][0].Indel); j++ {
					if dna.CompareSeqsIgnoreCase(alleles.Indel[i].Alt, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles.Indel[i].Ref, SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Alt) == 0 {
						SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Count = SampleMap[Location{loc.Chr, loc.Pos}][0].Indel[j].Count + alleles.Indel[i].Count
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
func ScoreVariants(input BatchSampleMap, sigThreshold float64, afThreshold float64) []*vcf.Vcf {

	fmt.Printf("#\n# Calling Variants\n")
	var VariantScores []*vcf.Vcf
	var progressMeter int
	var cA, cC, cG, cT, cIndel, dA, dC, dG, dT, dIndel int32
	var a, b []int32

	// Initialize a channel to send completed vcf structs through
	vcfChannel := make(chan *vcf.Vcf)
	var wg sync.WaitGroup

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
				a[i-1] = alleles[i].BaseA
				b[i-1] = alleles[0].BaseA - alleles[i].BaseA
			}
		case dna.C:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseC
				b[i-1] = alleles[0].BaseC - alleles[i].BaseC
			}
		case dna.G:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseG
				b[i-1] = alleles[0].BaseG - alleles[i].BaseG
			}
		case dna.T:
			for i = 1; i < len(alleles); i++ {
				a[i-1] = alleles[i].BaseT
				b[i-1] = alleles[0].BaseT - alleles[i].BaseT
			}
		default:
			continue
		}

		// Loop through samples and generate scores
		// For loop starts at index 1 because index zero is the background values
		for i = 1; i < len(alleles); i++ {

			// Retrieve Values for c
			cA = alleles[i].BaseA
			cC = alleles[i].BaseC
			cG = alleles[i].BaseG
			cT = alleles[i].BaseT

			// Retrieve Values for d
			dA = alleles[0].BaseA - alleles[i].BaseA
			dC = alleles[0].BaseC - alleles[i].BaseC
			dG = alleles[0].BaseG - alleles[i].BaseG
			dT = alleles[0].BaseT - alleles[i].BaseT

			// Generate Scores
			if alleles[i].Ref != dna.A {
				wg.Add(1)
				go score(&wg, vcfChannel, a[i-1], b[i-1], cA, dA, afThreshold, alleles[i], "A", 0, loc.Chr, loc.Pos, sigThreshold)
			}
			if alleles[i].Ref != dna.C {
				wg.Add(1)
				go score(&wg, vcfChannel, a[i-1], b[i-1], cC, dC, afThreshold, alleles[i], "C", 0, loc.Chr, loc.Pos, sigThreshold)
			}
			if alleles[i].Ref != dna.G {
				wg.Add(1)
				go score(&wg, vcfChannel, a[i-1], b[i-1], cG, dG, afThreshold, alleles[i], "G", 0, loc.Chr, loc.Pos, sigThreshold)
			}
			if alleles[i].Ref != dna.T {
				wg.Add(1)
				go score(&wg, vcfChannel, a[i-1], b[i-1], cT, dT, afThreshold, alleles[i], "T", 0, loc.Chr, loc.Pos, sigThreshold)
			}

			// Calculate p for each Indel
			for j = 0; j < len(alleles[i].Indel); j++ {
				cIndel = alleles[i].Indel[j].Count
				// Find Indel in the background Indel slice
				for l = 0; l < len(alleles[0].Indel); l++ {
					if dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Alt, alleles[0].Indel[l].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(alleles[i].Indel[j].Ref, alleles[0].Indel[l].Ref) == 0 {
						dIndel = alleles[0].Indel[l].Count - alleles[i].Indel[j].Count
						break
					}
				}
				wg.Add(1)
				go score(&wg, vcfChannel, a[i-1], b[i-1], cIndel, dIndel, afThreshold, alleles[i], "Indel", j, loc.Chr, loc.Pos, sigThreshold)
			}
		}
	}


	fmt.Println("# Waiting for Goroutines to finish")

	go func() {
		fmt.Println(runtime.NumGoroutine())
		wg.Wait()
		fmt.Println("# Goroutines finished")
		close(vcfChannel)
	}()

	for answer := range vcfChannel {
		VariantScores = append(VariantScores, answer)
	}

	return VariantScores
}

// Includes logic to exclude putative variants for which fishers exact test is unnecessary (e.g. alt allele count = 0) and exports as vcf.Vcf
func score(wg *sync.WaitGroup, vcfChannel chan *vcf.Vcf, a int32, b int32, c int32, d int32, afThreshold float64, inStruct *BatchAlleleCount, altbase string, indelslicepos int, chr string, pos int64, sigThreshold float64) {

	var p float64
	var answer *vcf.Vcf
	defer wg.Done()

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

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be exluded
	case float64(c)/float64(c+a) < afThreshold:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(int(a), int(b), int(c), int(d), true)
	}
	if p < sigThreshold {
		answer = &vcf.Vcf{
			Chr:     chr,
			Id:      ".",
			Qual:    p,
			Filter:  ".",
			Format:  "Sample:RefCount:AltCount:Cov",
			Unknown: fmt.Sprintf("%s:%d:%d:%d", inStruct.Sample, a, c, inStruct.Counts)}

		switch altbase {
		case "A", "C", "G", "T":
			answer.Pos = pos + 1
			answer.Info = "."
			answer.Ref = dna.BaseToString(inStruct.Ref)
			answer.Alt = altbase

		case "Indel":
			answer.Pos = pos
			answer.Info = "."
			answer.Ref = dna.BasesToString(inStruct.Indel[indelslicepos].Ref)
			answer.Alt = dna.BasesToString(inStruct.Indel[indelslicepos].Alt)
		}
		vcfChannel <- answer
	}
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

		Ukn := strings.Split(data[i].Unknown, ":")
		Info := strings.Split(data[i].Info, "=")

		if len(Info) < 2 {
			continue
		}

		// Eff Example Line: ANN=A|missense_variant|MODERATE|TGFBR3|TGFBR3|transcript|NM_003243.4|protein_coding|13/17|c.2050G>T|p.Asp684Tyr|2565/6467|2050/2556|684/851||
		Ann := strings.Split(Info[1], "|")
		Cov := strings.Split(Ukn[3], "\n")

		fmt.Fprintf(output, "%s,%s,%d,%s,%s,%.3v,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
			Ukn[0],       // Sample
			data[i].Chr,  // Chr
			data[i].Pos,  // Pos
			data[i].Ref,  // Ref
			data[i].Alt,  // Alt
			data[i].Qual, // pVal
			Ann[3],       // Gene
			Ann[9],       // DNA
			Ann[10],      // Protein
			Ukn[1],       // RefCount
			Ukn[2],       // AltCount
			Cov[0],       // Coverage
			Ann[1],       // Consequence
			Ann[2],       // Prediction
			Ann[6])       // Transcript
	}
}
