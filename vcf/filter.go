package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"strings"
)

func ASFilter(v *Vcf, AaAa []int16, BBbb []int16) bool {
	gt := genotypeHelper(v)
	if Unqualified(AaAa, gt) || Unqualified(BBbb, gt) {
		return false
	} else if Heterozygous(AaAa, gt) && UniqueHomozygous(BBbb, gt) {
		return true
	} else {
		return false
	}
}

func WrongFilter(v *Vcf, AaAa []int16, BBbb []int16) bool {
	gt := genotypeHelper(v)
	if Homozygous(AaAa, gt) || Heterozygous(BBbb, gt) {
		return true
	} else {
		return false
	}
}

func Homozygous(key []int16, gt []Haplotype) bool {
	for _, Aa := range key {
		if !IsHomozygous(gt[Aa]) {
			return false
		}
	}
	return true
}

func HetsOnly(v *Vcf, AaAa []int16, BBbb []int16) bool {
	gt := genotypeHelper(v)
	if Heterozygous(AaAa, gt) {
		return true
	} else {
		return false
	}
}

func WriteHeaderDev(file *os.File, header *VcfHeader) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
	}
	common.ExitIfError(err)
}

func SimpleASFilter(v *Vcf, axb int16, a int16, b int16) bool {
	gt := genotypeHelper(v)
	if LowQual(gt[axb]) || LowQual(gt[a]) || LowQual(gt[b]) {
		return false
	}
	if gt[axb].One != gt[axb].Two && (gt[a].One == gt[a].Two && (gt[a].One != gt[b].One || gt[a].Two != gt[b].Two)) && (gt[b].One == gt[b].Two) {
		return true
	} else {
		return false
	}
}
func ReadFilterList(filename string, dict map[string]int16) ([]string, []string) {
	file := fileio.Read(filename)
	if len(file) != 2 {
		log.Fatalf("Error: Sample sheet should be exactly 2 lines: hets one line 1, homo on line 2")
	}
	AaAa := strings.Split(file[0], ",")
	//log.Printf("read in %s\n", file[0])
	BBbb := strings.Split(file[1], ",")
	return AaAa, BBbb
	//log.Printf("read in %s\n", file[1])
	//log.Printf("mapped ids HET=%d, HOM=%d\n", MapNameToIndex(dict, AaAa), MapNameToIndex(dict, BBbb))
	//return MapNameToIndex(dict, AaAa), MapNameToIndex(dict, BBbb)
}

func MapNameToIndex(dict map[string]int16, list []string) []int16 {
	var answer []int16 = make([]int16, len(list))
	for i := 0; i < len(list); i++ {
		answer[i] = dict[list[i]]
	}
	return answer
}



/*
func FilterList(v *Vcf, hets []int16, homo []int16) bool {
	//Converts notes to slice of Haplotype structs
	gt := genotypeHelper(v)

	for i := 0 ; i <len(homo); i++ {

	}
}*/

/*
func MediumASFilter(v *Vcf, axb1 int16, axb2 int16, axb3 int16, a int16, b int16) bool {
	gt := genotypeHelper(v)
	if IsHeterozygous(gt[axb1]) || gt[axb3].One != gt[axb3].Two) && ((gt[a].One != gt[b].One && gt[a].Two != gt[b].Two)) && (gt[b].One == gt[b].Two) {
		return true
	} else {
		return false
	}
}*/

//all records must be Hets.
func StrictHetFilter(v *Vcf) bool {
	samples := genotypeHelper(v)
	for _, gVcf := range samples {
		if gVcf.One == gVcf.Two {
			return false
		}
	}
	return true
}

func MediumHetFilter(v *Vcf, index []int) bool {
	samples := genotypeHelper(v)
	for _, idx := range index {
		if samples[idx].One == samples[idx].Two {
			return false
		}
	}
	return true
}

func FilterAxtVcf(vcfs []*Vcf, fa []*fasta.Fasta) []*Vcf {
	split := VcfSplit(vcfs, fa)
	var answer []*Vcf
	var i, j int
	var ref []dna.Base
	var alt []dna.Base
	for i = 0; i < len(split); i++ {
		encountered := make(map[int64]bool)
		for j = 0; j < len(split[i]); j++ {
			if encountered[split[i][j].Pos] == true {
				//do not add
			} else {
				encountered[split[i][j].Pos] = true
				ref = dna.StringToBases(split[i][j].Ref)
				alt = dna.StringToBases(split[i][j].Alt)
				if dna.CountBaseInterval(ref, dna.N, 0, len(ref)) == 0 && dna.CountBaseInterval(alt, dna.N, 0, len(alt)) == 0 {
					answer = append(answer, split[i][j])
				}
			}
		}
	}
	Sort(answer)
	return answer
}

func FilterNs(vcfs []*Vcf) []*Vcf {
	var answer []*Vcf
	for i := 0; i < len(vcfs); i++ {
		if !strings.Contains(vcfs[i].Ref, "N") && !strings.Contains(vcfs[i].Alt, "N") {
			answer = append(answer, vcfs[i])
		}
	}
	return answer
}

func sameRecord(a *Vcf, b *Vcf) bool {
	if isEqual(a, b) {
		return true
	}
	if strings.Compare(a.Chr, b.Chr) == 0 && a.Pos == b.Pos {
		if strings.Compare(a.Ref, b.Ref) == 0 && strings.Compare(a.Alt, b.Alt) == 0 {
			return true
		}
	}
	return false
}

func mergeSimilarVcf(a *Vcf, b *Vcf) *Vcf {
	mergeRecord := &Vcf{Chr: a.Chr, Pos: a.Pos, Id: a.Id, Ref: "", Alt: "", Qual: a.Qual, Filter: "Merged:SNP:INDEL", Info: a.Info, Format: "SVTYPE=SNP", Notes: a.Notes}
	if len(a.Ref) < len(b.Ref) {
		mergeRecord.Ref += b.Ref
	} else {
		mergeRecord.Ref += a.Ref
	}
	if len(a.Alt) < len(b.Alt) {
		mergeRecord.Alt += b.Alt
	} else {
		mergeRecord.Alt += a.Alt
	}
	return mergeRecord
}
