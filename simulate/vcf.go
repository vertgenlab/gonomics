package simulate

import (
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/numbers"
)

//SimulateVCF generates simulated VCF data. This currently is a skeleton that can be expanded on with additional functions. Currently supports generating gVCF annotations from
//SimulateAFS in the popgen package. Random locations can be generated with simulateBed, and random mutations can be picked with Christi's simulate code.
//maybe we can combine these?
func SimulateVcf(alpha float64, n int, k int, outFile string) {
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	f := popgen.AFSStationarityClosure(alpha)
	var genotype []vcf.GenomeSample
	binHeights, sumHeights := numbers.InitializeFastRejectionSampler(0.001, 0.999, f, 1000)//set with hardcoded allele frequencies and 1000 bins for now.

	var current *vcf.Vcf
	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < k; i++ {
		genotype = popgen.SimulateGenotype(alpha, n, binHeights, sumHeights)
		//most fields are hardcoded but can be filled in later
		current = &vcf.Vcf{Chr: "chr1", Pos: int64(i+1), Id: ".", Ref: "A", Alt: "T", Qual: 100, Filter: ".", Info: ".", Format: ".", Notes: vcf.GenotypeToStringNew(genotype)}
		vcf.WriteVcf(out, current)
	}
}
