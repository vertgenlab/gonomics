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
	var curr *LinkedRead
	var done bool

	fwd := fileio.EasyOpen(fileOne)
	defer fwd.Close()

	rev := fileio.EasyOpen(fileTwo)
	defer rev.Close()

	for curr, done = NextSeq(fwd, rev); !done; curr, done = NextSeq(fwd, rev) {
		tenXg <- curr
	}
	close(tenXg)
}

func NextSeq(fwdReader *fileio.EasyReader, revReader *fileio.EasyReader) (*LinkedRead, bool) {
	curr := &PairedEnd{Fwd: nil, Rev: nil}
	var done bool
	curr, done = NextFastqPair(fwdReader, revReader)
	if curr == nil || done {
		return nil, true
	} else {
		//var tenXg *LinkedRead 
		//tenXg = 
		return FastqPairLinked(curr), false
	}
}

func FastqPairLinked(fqPair *PairedEnd) *LinkedRead {
	tenXG := LinkedRead{Fwd: nil, Rev: nil, Bx: GetBarcode(fqPair.Fwd, 0, 16), Umi: GetBarcode(fqPair.Fwd, 17,23)}
	tenXG.Fwd = TrimFastq(fqPair.Fwd, 24, len(fqPair.Fwd.Seq))
	tenXG.Rev = fqPair.Rev

	tenXG.Fwd.Name = fmt.Sprintf("%s %s-1", fqPair.Fwd.Name, dna.BasesToString(tenXG.Bx))
	tenXG.Rev.Name = fmt.Sprintf("%s %s-1", fqPair.Rev.Name, dna.BasesToString(tenXG.Bx))
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
/*
func readOneHandle(line1 string, line2 string, line3 string, line4 string) (*Fastq, []dna.Base, []dna.Base) {
	var curr Fastq

	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	curr = Fastq{Name: line1[1:len(line1)], Seq: dna.StringToBases(line2), Qual: []rune(line4)}
	bx := curr.Seq[:15]
	umi := curr.Seq[16:23]
	curr.Seq = curr.Seq[24:]
	return &curr, bx, umi
}

func LinkedROne(reader *fileio.EasyReader) (*Fastq, []dna.Base, []dna.Base, bool) {
	line, done := fileio.EasyNextLine(reader)
	line2, done2 := fileio.EasyNextLine(reader)
	line3, done3 := fileio.EasyNextLine(reader)
	line4, done4 := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	if done2 || done3 || done4 {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}
	return readOneHandle(line, line2, line3, line4), false
}

func LinkedReadToChan(readOne string, readTwo string, out chan<- *LinkedRead) {
	var readClouds *LinkedRead
	var done bool

	fileOne := fileio.EasyOpen(readOne)
	defer fileOne.Close()
	fileTwo := fileio.EasyOpen(readTwo)
	defer fileTwo.Close()

	for readClouds, done = NextLinked(fileOne, fileTwo); !done; readClouds, done = NextLinked(fileOne, fileTwo) {
		out <- readClouds
	}
	close(out)
}

func NextLinked(readerFwd *fileio.EasyReader, readerRev *fileio.EasyReader) (*LinkedRead, bool) {
	curr := LinkedRead{Fwd: nil, Rev: nil, Bx: nil, Umi: nil}
	fqOne, tenX, uniq, doneFwd := LinkedROne(readerFwd)
	fqTwo, doneRev := NextFastq(readerRev)

	curr.Fwd = fqOne
	curr.Rev = fqTwo
	curr.Bx = tenX
	curr.Umi = uniq
	return &curr, false
}*/

