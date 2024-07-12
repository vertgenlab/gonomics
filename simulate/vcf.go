package simulate

import (
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
)

// VcfToFile generates simulated VCF data.  The inputs are alpha (the selection parameter), the number of sites,
// the output filename, along with parameters for the bounding function for sampling.  Reasonable parameters
// choices for boundAlpha, boundBeta, and boundMultiplier are 0.001, 0.001, and 10000.
func VcfToFile(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64, setSeed int64) {
	out := fileio.EasyCreate(outFile)
	var current vcf.Vcf
	seed := rand.New(rand.NewSource(setSeed))

	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < numSites; i++ {
		current = SingleVcf(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, i+1, seed)
		vcf.WriteVcf(out, current)
	}

	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

// SingleVcf returns a single simulated Vcf record for a user-specified selection parameter alpha and genomic position.
// There also needs to be parameters for the bounding function, where alpha, beta, and multiplier parameters of 0.001, 0.001, and 10000 are good
// for most applications.
func SingleVcf(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int, seed *rand.Rand) vcf.Vcf {
	var genotype []vcf.Sample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, seed)
	//most fields are hardcoded but can be filled in later
	answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}
