package fastq

import(
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"fmt"

)
type LinkedRead struct {
	Fwd *Fastq
	Rev *Fastq
	Bx []dna.Base
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
	if curr == nil || done {
		return nil, true
	}
	return FastqPairLinked(curr), false
}

func FastqPairLinked(fqPair *PairedEnd) *LinkedRead {
	tenXG := LinkedRead{Fwd: nil, Rev: fqPair.Rev, Bx: GetBarcode(fqPair.Fwd, 0, 16), Umi: GetBarcode(fqPair.Fwd, 17,23)}
	tenXG.Fwd = TrimFastq(fqPair.Fwd, 24, len(fqPair.Fwd.Seq))
	bxTag := fmt.Sprintf("BX:%s", dna.BasesToString(tenXG.Bx))
	tenXG.Fwd.Name = fmt.Sprintf("%s_%s", fqPair.Fwd.Name, bxTag)
	tenXG.Rev.Name = fmt.Sprintf("%s_%s", fqPair.Rev.Name, bxTag)
	return &tenXG
}
//trimes fastq bases and quals given a start and end position zero base
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
//fastq files that have read1, read2 merged one after another
func InterLeaveFq(reader *fileio.EasyReader) (*PairedEnd, bool) {
	curr := PairedEnd{Fwd: nil, Rev: nil}
	fqOne, done1 := NextFastq(reader)
	fqTwo, done2 := NextFastq(reader)
	curr.Fwd = fqOne
	curr.Rev = fqTwo
	if done1 || done2 {
		return &curr, true
	}
	return &curr, false
}

func ReadInterLeaceLoop(filename string) []*PairedEnd {
	var curr *PairedEnd
	var done bool
	var answer []*PairedEnd
	er := fileio.EasyOpen(filename)
	for curr, done = InterLeaveFq(er); !done; curr, done = InterLeaveFq(er) {
		answer = append(answer, curr)
	}
	return answer
}

func fastqStats(fq *Fastq) {
	fmt.Printf("%s\t%d\n%s\n", fq.Name, len(fq.Seq), dna.BasesToString(fq.Seq))
}

func PrettyPrint(lr *LinkedRead) {
	fmt.Printf("Read\t%s\n10xG\t%s\nUmi\t%s\n\n", lr.Fwd.Name, dna.BasesToString(lr.Bx), dna.BasesToString(lr.Umi))
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
