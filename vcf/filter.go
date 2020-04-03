package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"strings"
)

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
