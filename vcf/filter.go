package vcf

import (
	"log"
	"strings"
)

func FilterAxtVcf(vcfs []*Vcf) []*Vcf {
	if len(vcfs) == 0 {
		log.Fatalf("Error: vcf file provided is empty...")
	}
	Sort(vcfs)
	var answer []*Vcf
	for i := 0; i < len(vcfs)-1; {
		if CompareVcf(vcfs[i], vcfs[i+1]) != 0 {
			answer = append(answer, vcfs[i])
			i++
		} else {
			//if can merge
			if vcfs[i].Pos == vcfs[i+1].Pos {
				if strings.Compare(vcfs[i].Ref, vcfs[i+1].Ref) == 0 || strings.Compare(vcfs[i].Alt, vcfs[i+1].Alt) == 0 {
					vcfs[i] = mergeSimilarVcf(vcfs[i], vcfs[i+1])
					answer = append(answer, vcfs[i])
					i += 2
				}
			}
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

/*
else if CompareName(curr, vcfs[i]) == 0 && curr.Pos+1 == vcfs[i].pos {
			if strings.Compare(curr.Format, "SVTYPE=SNP") == 0 && strings.Compare("SVTYPE=SNP", vcfs[i].Format) == 0 {
				for j = i; j < len(vcfs); j++ {
					curr.Ref+= vcfs[j].Ref
					curr.Alt+= vcfs[j].Alt
				}


			}
		}*/
