package simulate

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
)	

// VcfToFile generates simulated VCF data.  The inputs are alpha (the selection parameter), the number of sites,
// the output filename, along with parameters for the bounding function for sampling.  Reasonable parameters
// choices for boundAlpha, boundBeta, and boundMultiplier are 0.001, 0.001, and 10000.
func VcfToFile(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64, refFile string, hasRef bool) {
	out := fileio.EasyCreate(outFile)
	var current vcf.Vcf

	numGeneratedSites := 0
	var region bed.Bed
	var currKey int
	var foundInMap bool

	//for each segregating site, we make a vcf entry and write out
	if hasRef {
		generatedPos := make(map[int]bool)
		refFa := fasta.Read(refFile)
		bedRefFile := bed.UngappedRegionsAllFromFa(refFa)

		// changes bed region names to reference where they start/end at Ns
		regionOffset := mapSearchspaceToOffset(bedRefFile)
		faIndices := regionNameToFaIdx( refFa)

		var refBase dna.Base
		var regionNameSplit []string
		var regionNameStripped string
		totalWindows := CountWindows(bedRefFile, 1)
		var windowNumber int
		for numGeneratedSites < numSites {
			// generate a random position in the ungapped bed
			windowNumber = numbers.RandIntInRange(0, totalWindows) 
			region, _ = GenerateBedRegion(bedRefFile, windowNumber, 1)
			regionNameSplit = strings.Split(region.Name, "_") // split bed name to get fasta sequence original name
			regionNameStripped = regionNameSplit[0]
			
			currKey = regionOffset[regionNameStripped] + region.ChromStart
			// check if currKey not overlapping with previously generated positions
			if _, foundInMap = generatedPos[currKey]; !foundInMap {
				refBaseIndex := faIndices[regionNameStripped]
				refBase = refFa[refBaseIndex].Seq[region.ChromStart]
				current = SingleVcfWithRef(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, regionNameStripped, currKey+1, refBase)
				vcf.WriteVcf(out, current)
				numGeneratedSites++
				generatedPos[currKey] = true
			}
		}
	} else {
		for i := 0; i < numSites; i++ {
			current = SingleVcfRandom(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, i+1)
			vcf.WriteVcf(out, current)
		}
	}

	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

// regionNameToFaIdx maps sequence names in a fasta to their index
func regionNameToFaIdx(refFa []fasta.Fasta) map[string]int {
	faIndices := make(map[string]int)
	for idx, chrom := range refFa {
		faIndices[chrom.Name] = idx
	}
	return faIndices
}

// mapSearchspaceToOffset maps the start of each region in a bed to the overall position relative to the first region appearing in the bed
func mapSearchspaceToOffset(searchSpace []bed.Bed) map[string]int {
	regionOffset := make(map[string]int)
	prevRegionEnd := 0
	for _, region := range searchSpace {
		regionOffset[region.Name] = prevRegionEnd
		prevRegionEnd = prevRegionEnd + region.ChromEnd
	}
	return regionOffset
}

// SingleVcfRandom returns a single simulated Vcf record for a user-specified selection parameter alpha and genomic position.
// There also needs to be parameters for the bounding function, where alpha, beta, and multiplier parameters of 0.001, 0.001, and 10000 are good
// for most applications.
func SingleVcfRandom(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int) vcf.Vcf {
	var genotype []vcf.Sample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)

	// ref := ChooseRandomBase(float64(0.42))
	// answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: dna.BaseToString(ref), Alt: []string{dna.BaseToString(changeBase(ref))}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	
	// keeping hardcoded for old tests
	answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}

	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}

// SingleVcfWithRef returns a single simulated Vcf record for a user-specified selection parameter alpha, genomic position and reference Fasta file.
// There also needs to be parameters for the bounding function, where alpha, beta, and multiplier parameters of 0.001, 0.001, and 10000 are good
// for most applications.
func SingleVcfWithRef(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, chrom string, pos int, refBase dna.Base) vcf.Vcf {
	var genotype []vcf.Sample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)

	answer = vcf.Vcf{Chr: chrom, Pos: pos, Id: ".", Ref: dna.BaseToString(refBase), Alt: []string{dna.BaseToString(changeBase(refBase))}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}
