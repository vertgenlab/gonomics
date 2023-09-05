package vcf

import (
	"log"
	//DEBUG: "fmt".
)

// invertSample inverts the ancestral/derived state for each allele in a Sample. Only works for biallelic positions, throws an error if an allele state is greater than 1.
// TODO: this is currently only supported for biallelic bases.
func invertSample(g Sample) Sample {
	for i := range g.Alleles {
		switch g.Alleles[i] {
		case 0:
			g.Alleles[i] = 1
		case 1:
			g.Alleles[i] = 0
		default:
			log.Fatal("Error in InvertGenomeSample: bases must be biallelic to be inverted.")
		}
	}
	return g
}

// invertAlleles inverts the Genotype for each entry in a slice of Sample structs.
func invertAlleles(g []Sample) []Sample {
	for i := 0; i < len(g); i++ {
		g[i] = invertSample(g[i])
	}
	return g
}

// InvertVcf inverts the reference and alt variants in a Vcf record. Currently does not update other fields, but this functionality may be added.
func InvertVcf(v Vcf) Vcf {
	if len(v.Alt) > 1 {
		log.Fatalf("InvertVCF is not currently supported for multiallelic bases.")
	}
	//DEBUG: fmt.Printf("v.Ref before inversion: %s.\n", v.Ref)
	v.Ref, v.Alt[0] = v.Alt[0], v.Ref
	v.Samples = invertAlleles(v.Samples)
	//DEBUG: fmt.Printf("v.Ref after inversion: %s.\n", v.Ref)
	return v
}
