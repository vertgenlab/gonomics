package simulate

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
)

//SimulateVCF generates simulated VCF data. This currently is a skeleton that can be expanded on with additional functions. Currently supports generating gVCF annotations from
//SimulateAFS in the popgen package. Random locations can be generated with simulateBed, and random mutations can be picked with Christi's simulate code.
//maybe we can combine these?
func SimulateVcf(alpha float64, n int, k int, outFile string) {
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	var genotype []vcf.GenomeSample

	var current *vcf.Vcf
	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < k; i++ {
		genotype = popgen.SimulateGenotype(alpha, n)
		//most fields are hardcoded but can be filled in later
		current = &vcf.Vcf{Chr: "chr1", Pos: i+1, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}

		vcf.WriteVcf(out, current)
	}
}
