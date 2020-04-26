package vcf

import(
	"strings"
	"fmt"
)

func ShuffleVcfColumn(input *Vcf, het[]int16, hom []int16) *Vcf {
	columnData := make([]string, 0,len(het)+len(hom))
	words := strings.Split(input.Notes, "\t")
	var i int
	for i = 0; i < len(het); i++ {
		columnData = append(columnData, words[het[i]])
	}
	for i = 0; i < len(hom);i++ {
		columnData = append(columnData, words[hom[i]])
	}
	input.Notes = strings.Join(columnData, "\t")
	return input
}

func ViewGenotypeVcf(v *Vcf) {
	gVcf := vcfToGenotype(v)
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, haplotypesToString(gVcf))
}

func prettyShuffle(v *Vcf, het[]int16, hom []int16) {
	gVcf := vcfToGenotype(ShuffleVcfColumn(v, het, hom))
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, haplotypesToString(gVcf))	
}

func vcfPrettyPrint(v *Vcf) {
		fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, v.Notes)
}