package vcf

import (
	"fmt"
	"strings"
)

func ShuffleVcfColumn(input *Vcf, het []int16, hom []int16) *Vcf {
	columnData := make([]string, 0, len(het)+len(hom))
	words := strings.Split(input.Notes, "\t")
	var i int
	for i = 0; i < len(het); i++ {
		columnData = append(columnData, words[het[i]])
	}
	for i = 0; i < len(hom); i++ {
		columnData = append(columnData, words[hom[i]])
	}
	input.Notes = strings.Join(columnData, "\t")
	return input
}

func ViewGenotypeVcf(v *Vcf) {
	gVcf := VcfToGenotype(v)
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, haplotypesToString(gVcf))
}

func PrettyShuffle(v *Vcf, het []int16, hom []int16) {
	gVcf := VcfToGenotype(ShuffleVcfColumn(v, het, hom))
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, haplotypesToString(gVcf))
}

func vcfPrettyPrint(v *Vcf) {
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, v.Notes)
}
