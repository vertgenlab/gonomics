package gtf

import (
	"strings"
)

func AllAreEqual(a map[string]*Gene, b map[string]*Gene) bool {
	if len(a) != len(b) {
		return false
	}
	for k, _ := range b {
		if _, ok := a[k]; !ok {
			return false
		} else {
			if !EqualGene(a[k], b[k]) {
				return false
			}
		}
	}
	return true
}

func EqualGene(a *Gene, b *Gene) bool {
	if strings.Compare(a.GeneID, b.GeneID) != 0 {
		return false
	}
	if strings.Compare(a.GeneName, b.GeneName) != 0 {
		return false
	}
	for i := 0; i < len(a.Transcripts); i++ {
		if !EqualTranscript(a.Transcripts[i], b.Transcripts[i]) {
			return false
		}
	}
	return true
}

func EqualTranscript(a *Transcript, b *Transcript) bool {
	if strings.Compare(a.Chr, b.Chr) != 0 {
		return false
	}
	if strings.Compare(a.Source, b.Source) != 0 {
		return false
	}
	if a.Start != b.Start || a.End != b.End {
		return false
	}
	if a.Score != b.Score {
		return false
	}
	if a.Strand != b.Strand {
		return false
	}
	if strings.Compare(a.TranscriptID, b.TranscriptID) != 0 {
		return false
	}

	for i := 0; i < len(a.Exons); i++ {
		if !EqualExon(a.Exons[i], b.Exons[i]) {
			return false
		}
	}
	return true
}

func EqualExon(a *Exon, b *Exon) bool {
	if strings.Compare(a.ExonID, b.ExonID) != 0 {
		return false
	}
	if strings.Compare(a.ExonNumber, b.ExonNumber) != 0 {
		return false
	}
	if a.Start != b.Start || a.End != b.End {
		return false
	}
	if a.Score != b.Score {
		return false
	}
	if a.FiveUtr != nil { //if a has a UTR we compare to b
		if b.FiveUtr == nil {
			return false
		}
		if !EqualFiveUtr(a.FiveUtr, b.FiveUtr) {
			return false
		}
	} else if b.FiveUtr != nil {
		return false
	}
	if a.Cds != nil {
		if b.Cds == nil {
			return false
		}
		if !EqualCds(a.Cds, b.Cds) {
			return false
		}
	} else if b.Cds != nil {
		return false
	}
	if a.ThreeUtr != nil {
		if b.ThreeUtr == nil {
			return false
		}
		if !EqualThreeUtr(a.ThreeUtr, b.ThreeUtr) {
			return false
		}
	} else if b.ThreeUtr != nil {
		return false
	}
	return true
}

func EqualFiveUtr(a *FiveUTR, b *FiveUTR) bool {
	return a.Start == b.Start && a.End == b.End && a.Score == b.Score
}

func EqualCds(a *CDS, b *CDS) bool {
	return a.Start == b.Start && a.End == b.End && a.Score == b.Score && a.Frame == b.Frame
}

func EqualThreeUtr(a *ThreeUTR, b *ThreeUTR) bool {
	return a.Start == b.Start && a.End == b.End && a.Score == b.Score
}
