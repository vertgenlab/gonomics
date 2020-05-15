package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
)

type LinkedRead struct {
	Fwd *Fastq
	Rev *Fastq
	Bx  []dna.Base
	Umi []dna.Base
}

func ReadToChanLinked(fileOne string, fileTwo string, tenXg chan<- *LinkedRead) {
	fwd, rev := fileio.EasyOpen(fileOne), fileio.EasyOpen(fileTwo)
	defer fwd.Close()
	defer rev.Close()

	for curr, done := NextSeq(fwd, rev); !done; curr, done = NextSeq(fwd, rev) {
		tenXg <- curr
	}
	close(tenXg)
}

func NextSeq(fwdReader *fileio.EasyReader, revReader *fileio.EasyReader) (*LinkedRead, bool) {
	curr, done := NextFastqPair(fwdReader, revReader)
	if done {
		return nil, true
	}
	return FastqPairLinked(curr), false
}

//10x linked-read library construct are modifications only on forward or read one, read two remains unchanged:
//1) first 16 bases: 10x barcode labeled either Bx or Rx tags
//2) Next 6 bases is a 6 base Umi
//3) finally adaptors and genomic sequence
//TODO: consider trimming off adapter sequences too
func FastqPairLinked(fqPair *PairedEnd) *LinkedRead {
	tenXG := LinkedRead{Fwd: nil, Rev: fqPair.Rev, Bx: GetBarcode(fqPair.Fwd, 0, 16), Umi: GetBarcode(fqPair.Fwd, 17, 23)}
	tenXG.Fwd = TrimFastq(fqPair.Fwd, 24, len(fqPair.Fwd.Seq))
	bxTag := fmt.Sprintf("BX:%s", dna.BasesToString(tenXG.Bx))
	tenXG.Fwd.Name = fmt.Sprintf("%s_%s", fqPair.Fwd.Name, bxTag)
	tenXG.Rev.Name = fmt.Sprintf("%s_%s", fqPair.Rev.Name, bxTag)
	return &tenXG
}

//trims fastq bases and quals given a start and end position zero base
func TrimFastq(fq *Fastq, start int, end int) *Fastq {
	fq.Seq = fq.Seq[start:end]
	fq.Qual = fq.Qual[start:end]
	return fq
}

//copies a subset of the given fastq record, at start and and end sites, zero base
//could be a UMI or 10x barcode
func GetBarcode(fq *Fastq, start int, end int) []dna.Base {
	answer := make([]dna.Base, end-start)
	copy(answer, fq.Seq[start:end])
	return answer
}

/*
//fastq files that have read1, read2 merged one after another
func InterLeaveFq(reader *fileio.EasyReader) (*PairedEnd, bool) {
	fqOne, done1 := NextFastq(reader)
	fqTwo, done2 := NextFastq(reader)
	if done1 || done2 {
		return nil, true
	}
	curr := PairedEnd{Fwd: nil, Rev: nil}
	curr.Fwd = fqOne
	curr.Rev = fqTwo
	return &curr, false
}

//Paired reads that are contained in one file, in alt. order (i.e. 1,2,1,2,1)
func ReadInterLeaceLoop(filename string) []*PairedEnd {
	var curr *PairedEnd
	var done bool
	var answer []*PairedEnd
	er := fileio.EasyOpen(filename)
	for curr, done = InterLeaveFq(er); !done; curr, done = InterLeaveFq(er) {
		answer = append(answer, curr)
	}
	return answer
}*/

func fastqStats(fq *Fastq) string {
	return fmt.Sprintf("%s\t%d\n%s\n", fq.Name, len(fq.Seq), dna.BasesToString(fq.Seq))
}

func PrettyPrint(lr *LinkedRead) string {
	return fmt.Sprintf("Read\t%s\n10xG\t%s\nUmi\t%s\n\n", lr.Fwd.Name, dna.BasesToString(lr.Bx), dna.BasesToString(lr.Umi))
}

//TODO: add a more stringent test to make sure we are trimming the reads correctly
func isLinkedRead(lr *LinkedRead) bool {
	if len(lr.Bx) != 16 {
		return false
	}
	if len(lr.Umi) != 6 {
		return false
	}
	return true
}
