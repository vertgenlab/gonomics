package simulate

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
)

//VcfToFile generates simulated VCF data. This currently is a skeleton that can be expanded on with additional functions. Currently supports generating gVCF annotations from
//SimulateAFS in the popgen package. Random locations can be generated with simulateBed, and random mutations can be picked with Christi's simulate code.
//maybe we can combine these?
func VcfToFile(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64) {
	out := fileio.EasyCreate(outFile)
	var current vcf.Vcf

	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < numSites; i++ {
		current = SingleVcf(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, i+1)
		vcf.WriteVcf(out, current)
	}

	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

//SingleVcf returns a single simulated Vcf from a user-specified selection parameter alpha.
func SingleVcf(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int) vcf.Vcf {
	var genotype []vcf.GenomeSample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)
	//most fields are hardcoded but can be filled in later
	answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}
