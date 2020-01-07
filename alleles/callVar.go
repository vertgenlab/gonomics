package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"io"
	"os"
)

type VarScore struct {
	Sample	string
	Ref		dna.Base
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
	//TODO: remove background values once troubleshooting is done
	BkgdA	int32
	BkgdC	int32
	BkgdG	int32
	BkgdT	int32
	BkgdIns	int32
	BkgdDel	int32
	pA 		float64
	pC 		float64
	pG 		float64
	pT 		float64
	pIns 	float64
	pDel 	float64
}

// Map structure: map[Chromosome]map[Position][Sample]*VarScore
type VarMap map[string]map[int64][]*VarScore

func ScoreVariants (input BatchSampleMap, sigThreshold float64) VarMap {

	fmt.Printf("#\n# Calling Variants\n")
	var current *VarScore
	VariantScores := make(map[string]map[int64][]*VarScore)
	var progressMeter int

	for chrName, chr := range input {
		for pos, alleles := range chr {

			if progressMeter % 50000 == 0 {
				fmt.Printf("# Processed %d Positions\n", progressMeter)
			}
			progressMeter++

			// Must be at least 2 samples in the slice to generate a batch p value.
			// A single sample would be len = 2 because the background values are always element 0
			if len(alleles) <= 2 {
				continue
			}

			//if the chromosome has already been added to the matrix, move along
			_, ok := VariantScores[chrName]

			//if the chromosome is NOT in the matrix, initialize
			if ! ok {
				VariantScores[chrName] = make(map[int64][]*VarScore)
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
			var i int
			var a []int32 = make([]int32, len(alleles)-1)
			var b []int32 = make([]int32, len(alleles)-1)
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
				cA := alleles[i].BaseA
				cC := alleles[i].BaseC
				cG := alleles[i].BaseG
				cT := alleles[i].BaseT
				cIns := alleles[i].Ins
				cDel := alleles[i].Del

				// Retrieve Values for d
				dA := alleles[0].BaseA - alleles[i].BaseA
				dC := alleles[0].BaseC - alleles[i].BaseC
				dG := alleles[0].BaseG - alleles[i].BaseG
				dT := alleles[0].BaseT - alleles[i].BaseT
				dIns := alleles[0].Ins - alleles[i].Ins
				dDel := alleles[0].Del - alleles[i].Del

				// Generate Scores
				pA := Score(a[i-1], b[i-1], cA, dA)
				pC := Score(a[i-1], b[i-1], cC, dC)
				pG := Score(a[i-1], b[i-1], cG, dG)
				pT := Score(a[i-1], b[i-1], cT, dT)
				pIns := Score(a[i-1], b[i-1], cIns, dIns)
				pDel := Score(a[i-1], b[i-1], cDel, dDel)

				// If any p value is below the significance threshold then append to the map
				pMin := MinScore(pA, pC, pG, pT, pIns, pDel)
				if pMin <= sigThreshold {
					current = &VarScore{
						Sample:		alleles[i].Sample,
						Ref:		alleles[i].Ref,
						BaseA:		alleles[i].BaseA,
						BaseC:		alleles[i].BaseC,
						BaseG:		alleles[i].BaseG,
						BaseT:		alleles[i].BaseT,
						Ins:		alleles[i].Ins,
						Del:		alleles[i].Del,
						//TODO: remove background values once troubleshooting is done
						BkgdA:		alleles[0].BaseA - alleles[i].BaseA,
						BkgdC:		alleles[0].BaseC - alleles[i].BaseC,
						BkgdG:		alleles[0].BaseG - alleles[i].BaseG,
						BkgdT:		alleles[0].BaseT - alleles[i].BaseT,
						BkgdIns:	alleles[0].Ins - alleles[i].Ins,
						BkgdDel:	alleles[0].Del - alleles[i].Del,
						pA: 		pA,
						pC: 		pC,
						pG: 		pG,
						pT: 		pT,
						pIns: 		pIns,
						pDel:		pDel}
					VariantScores[chrName][pos] = append(VariantScores[chrName][pos], current)
				}
			}
		}
	}

	return VariantScores
}


func Score (a int32, b int32, c int32, d int32) float64 {

	var p float64

	switch {
	// If alternate allele is zero then there is no variant and score is 1
	case c == 0:
		return 1

	// If a = b and c = d then it is testing itself and should return 1
	case a == b && c == d:
		return 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(int(a), int(b), int(c), int(d), true)
	}
	return p
}


func MinScore (pA float64, pC float64, pG float64, pT float64, pIns float64, pDel float64) float64 {
	min := pA
	if min > pC {min = pC}
	if min > pG {min = pG}
	if min > pT {min = pT}
	if min > pIns {min = pIns}
	if min > pDel {min = pDel}
	return min
}


func WriteVarMap (input VarMap, output string) {

	var outFile *os.File

	if output != "stdout" {
		fmt.Printf("#Creating Output File\n")
		outFile, _ = os.Create(output)
		defer outFile.Close()
		//TODO: remove background values once troubleshooting is done
		io.WriteString(outFile, "Sample\tChr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\tBkgdA\tBkgdC\tBkgdG\tBkgdT\tBkgdIns\tBkgdDel\tpA\tpC\tpG\tpT\tpIns\tpDel\n")
	} else {
		fmt.Printf("#No Output Specified. Printing to STDOUT.\n")
		fmt.Printf("Sample\tChr\tPos\tRef\tA\tC\tG\tT\tIns\tDel\tBkgdA\tBkgdC\tBkgdG\tBkgdT\tBkgdIns\tBkgdDel\tpA\tpC\tpG\tpT\tpIns\tpDel\n")
	}

	var base string
	var progressMeter int
	var i int

	for chrName, chr := range input {
		for pos, alleles := range chr {
			for i = 0; i < len(alleles); i++ {

				switch alleles[i].Ref {
				case dna.A:
					base = "A"
				case dna.C:
					base = "C"
				case dna.G:
					base = "G"
				case dna.T:
					base = "T"
				case dna.N:
					base = "N"
				case dna.Gap:
					base = "Gap"
				case dna.Dot:
					base = "Dot"
				default:
					base = "NA"
				}

				if output != "stdout" {

					if progressMeter % 500000 == 0 {
						fmt.Printf("Wrote %d Positions\n", progressMeter)
					}
					progressMeter++

					// Write to file
					fmt.Fprintf(outFile,
						"%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n",
						alleles[i].Sample,
						chrName,
						pos + 1,
						base,
						alleles[i].BaseA,
						alleles[i].BaseC,
						alleles[i].BaseG,
						alleles[i].BaseT,
						alleles[i].Ins,
						alleles[i].Del,
						//TODO: remove background values once troubleshooting is done
						alleles[i].BkgdA,
						alleles[i].BkgdC,
						alleles[i].BkgdG,
						alleles[i].BkgdT,
						alleles[i].BkgdIns,
						alleles[i].BkgdDel,
						alleles[i].pA,
						alleles[i].pC,
						alleles[i].pG,
						alleles[i].pT,
						alleles[i].pIns,
						alleles[i].pDel)
				} else {
					// Write to stdout
					fmt.Printf(
						"%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
						alleles[i].Sample,
						chrName,
						pos + 1,
						base,
						alleles[i].BaseA,
						alleles[i].BaseC,
						alleles[i].BaseG,
						alleles[i].BaseT,
						alleles[i].Ins,
						alleles[i].Del,
						//TODO: remove background values once throubleshooting is done
						alleles[i].BkgdA,
						alleles[i].BkgdC,
						alleles[i].BkgdG,
						alleles[i].BkgdT,
						alleles[i].BkgdIns,
						alleles[i].BkgdDel,
						alleles[i].pA,
						alleles[i].pC,
						alleles[i].pG,
						alleles[i].pT,
						alleles[i].pIns,
						alleles[i].pDel)
				}
			}
		}
	}
}