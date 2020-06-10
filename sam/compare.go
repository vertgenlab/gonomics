package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"sort"
	"strings"
)

func compareName(alpha *SamAln, beta *SamAln) int {
	return strings.Compare(alpha.RName, beta.RName)
}

func SortByName(alignments []*SamAln) {
	sort.Slice(alignments, func(i, j int) bool { return compareName(alignments[i], alignments[j]) == -1 })
}

func SortByCoord(samRecords []*SamAln) {
	sort.Slice(samRecords, func(i, j int) bool { return Compare(samRecords[i], samRecords[j]) == -1 })
}

func Compare(a *SamAln, b *SamAln) int {
	chrName := strings.Compare(a.RName, b.RName)
	if chrName != 0 {
		return chrName
	}
	if a.Pos < b.Pos {
		return -1
	}
	if a.Pos > b.Pos {
		return 1
	}
	return 0
}

func IsEqual(a *SamAln, b *SamAln) bool {
	if strings.Compare(a.QName, b.QName) != 0 {
		log.Printf("%s != %s", a.QName, b.QName)
		return false
	}
	if a.Flag != b.Flag {
		log.Printf("%d != %d", a.Flag, b.Flag)
		return false
	}
	if strings.Compare(a.RName, b.RName) != 0 {
		log.Printf("%s != %s", a.RName, b.RName)
		return false
	}
	if a.Pos != b.Pos {
		log.Printf("%d != %d", a.Pos, b.Pos)
		return false
	}
	if a.MapQ != b.MapQ {
		log.Printf("%d != %d", a.MapQ, b.MapQ)
		return false
	}
	if strings.Compare(cigar.ToString(a.Cigar), cigar.ToString(b.Cigar)) != 0 {
		log.Printf("%s != %s", cigar.ToString(a.Cigar), cigar.ToString(b.Cigar))
		return false
	}
	if strings.Compare(a.RNext, b.RNext) != 0 {
		log.Printf("%s != %s", a.RNext, b.RNext)
		return false
	}
	if a.PNext != b.PNext {
		log.Printf("%d != %d", a.PNext, b.PNext)
		return false
	}
	if a.TLen != b.TLen {
		log.Printf("%d != %d", a.TLen, b.TLen)
		return false
	}
	if dna.CompareSeqsIgnoreCase(a.Seq, b.Seq) != 0 {
		log.Printf("%s != %s", dna.BasesToString(a.Seq), dna.BasesToString(b.Seq))
		return false
	}
	if strings.Compare(a.Qual, b.Qual) != 0 {
		log.Printf("%s != %s", a.Qual, b.Qual)
		return false
	}
	if strings.Compare(a.Extra, b.Extra) != 0 {
		log.Printf("%s != %s", a.Extra, b.Extra)
		return false
	}
	return true
}

func SplitSamByChr(samRecords []*SamAln) map[string][]*SamAln {
	SortByCoord(samRecords)
	genome := make(map[string][]*SamAln)
	for _, alignedRead := range samRecords {
		genome[alignedRead.RName] = append(genome[alignedRead.RName], alignedRead)
	}
	return genome
}
