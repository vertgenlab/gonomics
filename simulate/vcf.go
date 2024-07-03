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
// choices for boundAlpha, boundBeta, and boundMultiplier are 0.001, 0.001, and 10000.
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

// SingleVcf returns a single simulated Vcf record for a user-specified selection parameter alpha and genomic position.
// There also needs to be parameters for the bounding function, where alpha, beta, and multiplier parameters of 0.001, 0.001, and 10000 are good
// for most applications.
func SingleVcf(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int) vcf.Vcf {
	var genotype []vcf.Sample
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

// VcfToFileWithFasta generates simulated VCF data.  The inputs are alpha (the selection parameter), the number of sites,
// the output filename, along with parameters for the bounding function for sampling.  Reasonable parameters
// choices for boundAlpha, boundBeta, and boundMultiplier are 0.001, 0.001, and 10000.
func VcfToFileWithFasta(alpha float64, numAlleles int, numSites int, outFile string, boundAlpha float64, boundBeta float64, boundMultiplier float64, refFile fasta.Fasta, hasRef bool) {
	out := fileio.EasyCreate(outFile)
	var current vcf.Vcf

	//for each segregating site, we make a vcf entry and write out
	if hasRef {
		for i := 0; i < numSites; i++ {
			current = SingleVcfWithRef(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier, i+1, refBase, hasRef)
			vcf.WriteVcf(out, current)
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

func SingleVcfWithRef(alpha float64, numAlleles int, boundAlpha float64, boundBeta float64, boundMultiplier float64, pos int, refBase dna.base) vcf.Vcf {
	var genotype []vcf.Sample
	var divergent bool
	var answer vcf.Vcf
	genotype, divergent = popgen.SimulateGenotype(alpha, numAlleles, boundAlpha, boundBeta, boundMultiplier)
	//most fields are hardcoded but can be filled in later
	// already exists: ability to simulte llele frequencies based on allele distribution (boundalpha, boundbeta, boundmultiplier)
	// using jukes-cantor (equal probability mutate to any other)
	answer = vcf.Vcf{Chr: "chr1", Pos: pos, Id: ".", Ref: dna.BaseToString(refBase), Alt: []string{dna.BaseToString(changeBase(refBase))}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}
	if divergent {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Alt[0]))
	} else {
		answer = vcf.AppendAncestor(answer, dna.StringToBases(answer.Ref))
	}
	return answer
}
