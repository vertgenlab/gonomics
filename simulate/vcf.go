package simulate

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/"
)

// VcfToFile generates simulated VCF data.  The inputs are alpha (the selection parameter), the number of sites,
// the output filename, along with parameters for the bounding function for sampling.  Reasonable parameters
// choices for boundAlpha, boundBeta, and boundMultiplier are 0.001, 0.001, and 10000. Specify a reference file to base the simulated VCF off on.
func VcfToFile(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64, refFile fasta.Fasta, hasRef bool) {
	out := fileio.EasyCreate(outFile)
	var current vcf.Vcf
	
	numGeneratedSites := 0
	var region bed.bed
	var currKey int
	var foundInMap bool
	//for each segregating site, we make a vcf entry and write out
	// cmt: use faformat to make nogapbed from input fasta ungappdregionsallfromfa in bed/info.go
	if hasRef {
		var positions map[int]bool
		refFa = fasta.Read(refFile)
		convertedRefFile := UngappedRegionsFromFa(refFa)
		var length, totalWindows int

		for numGeneratedSites < numSites {
			// generate random position
			length, totalWindows = CountWindows(convertedRefFile, 1)
			region = GenerateBedRegion(convertedRefFile, totalWindows, 1)
			
			// cmt: region.ChromStart will map to same place on different chroms
			// helperfunction to map chrom pos to cumulative pos
			// end cmt
			currKey = helperFunction(region)
			// check if not overlapping with previously generated positions

			if _, foundInMap = m[currKey]; !foundInMap {
				current = SingleVcfWithRef(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, region.chrom, region.ChromStart+1, refBase, hasRef)
				vcf.WriteVcf(out, current)
				numGeneratedSites++
				m[currKey] = true
			}
		}
	} else {
		for i := 0; i < numSites; i++ {
			current = SingleVcf(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, i+1)
			vcf.WriteVcf(out, current)
		}
	}

	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

// SingleVcf returns a single simulated Vcf record for a user-specified selection parameter alpha and genomic position.
// There also needs to be parameters for the bounding function, where alpha, beta, and multiplier parameters of 0.001, 0.001, and 10000 are good
// for most applications.
func SingleVcf(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int) vcf.Vcf {
	var genotype []vcf.Sample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)

	answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}

// TODO: comment this function
func SingleVcfWithRef(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, chrom string, pos int, refBase dna.base) vcf.Vcf {
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
