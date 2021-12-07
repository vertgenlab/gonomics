package simulate

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
)

//SimulateVCF generates simulated VCF data. This currently is a skeleton that can be expanded on with additional functions. Currently supports generating gVCF annotations from
//SimulateAFS in the popgen package. Random locations can be generated with simulateBed, and random mutations can be picked with Christi's simulate code.
//maybe we can combine these?
func Vcf(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64) {
	out := fileio.EasyCreate(outFile)
	var genotype []vcf.Sample

	var current vcf.Vcf
	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < numSites; i++ {
		genotype = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)
		//most fields are hardcoded but can be filled in later
		current = vcf.Vcf{Chr: "chr1", Pos: i + 1, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
		vcf.WriteVcf(out, current)
	}
	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}
