package vcf

import (
	"fmt"
	"testing"

	"log"
)

func TestGvcf(t *testing.T) {
	//aXb := Read("litcMataF1.het.vcf")
	//m := GenotypeToMap(aXb)
	vcfPipe := make(chan *Vcf)

	dict := HeaderToMaps("testdata/smallLITCxMATAGenotypeGVCFs.vcf")
	log.Printf("Found %d contigs in this header\n", len(dict.Fa))
	fmt.Printf("Sample index is: %v\n", dict.HapIdx)
	go ReadToChan("testdata/smallLITCxMATAGenotypeGVCFs.vcf", vcfPipe)
	for record := range vcfPipe {
		fmt.Printf("%s\n", genotypeToString(vcfToGenotype(record)))
	}


}
