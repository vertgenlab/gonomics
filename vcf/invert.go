package vcf

import (
	"log"
	//DEBUG: "fmt"
)

//InvertGenomeSample inverts the ancestral/derived state for each allele in a Sample. Only works for biallelic positions, throws an error if an allele state is greater than 1.
//TODO: this is currently only supported for biallelic bases.
func InvertGenomeSample(g Sample) Sample {
	if g.AlleleOne == 0 {
		g.AlleleOne = 1
	} else if g.AlleleOne == 1 {
		g.AlleleOne = 0
	} else {
		log.Fatalf("Error in InvertGenomeSample: bases must be biallele to be inverted.")
	}
	if g.AlleleTwo == 0 {
		g.AlleleTwo = 1
	} else if g.AlleleTwo == 1 {
		g.AlleleTwo = 0
	} else {
		log.Fatalf("Error in InvertGenomeSample: bases must be biallele to be inverted.")
	}

	return g
}

//InvertAlleles inverts the Genotype for each entry in a slice of Sample structs.
func InvertAlleles(g []Sample) []Sample {
	for i := 0; i < len(g); i++ {
		g[i] = InvertGenomeSample(g[i])
	}
	return g
}

//InvertVcf inverts the reference and alt variants in a Vcf record. Currently does not update other fields, but this functionality may be added.
func InvertVcf(v Vcf) Vcf {
	if len(v.Alt) > 1 {
		log.Fatalf("InvertVCF is not currently supported for polyallelic bases.")
	}
	//DEBUG: fmt.Printf("v.Ref before inversion: %s.\n", v.Ref)
	v.Ref, v.Alt[0] = v.Alt[0], v.Ref
	v.Samples = InvertAlleles(v.Samples)
	//DEBUG: fmt.Printf("v.Ref after inversion: %s.\n", v.Ref)
	return v
}
