package fastq

import (
	"github.com/vertgenlab/gonomics/fileio"
	//"log"
	//"strings"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"io"
)

type PairedEnd struct {
	Fwd *Fastq
	Rev *Fastq
}

func ReadPairs(readOne string, readTwo string) []*PairedEnd {
	file1 := fileio.EasyOpen(readOne)
	defer file1.Close()

	file2 := fileio.EasyOpen(readTwo)
	defer file2.Close()

	answer := ReadFastqsPairs(file1, file2)
	return answer
}

func PairEndToChan(readOne string, readTwo string, output chan<- *PairedEnd) {
	var curr *PairedEnd
	var done bool

	fileOne := fileio.EasyOpen(readOne)
	defer fileOne.Close()
	fileTwo := fileio.EasyOpen(readTwo)
	defer fileTwo.Close()

	for curr, done = NextFastqPair(fileOne, fileTwo); !done; curr, done = NextFastqPair(fileOne, fileTwo) {
		output <- curr
	}
	close(output)
}

//TODO: rewrite logic to catch error
func NextFastqPair(reader1 *fileio.EasyReader, reader2 *fileio.EasyReader) (*PairedEnd, bool) {
	curr := PairedEnd{Fwd: nil, Rev: nil}
	fqOne, done1 := NextFastq(reader1)
	fqTwo, done2 := NextFastq(reader2)
	//name1 := strings.Split(fqOne.Name, " ")
	//name2 := strings.Split(fqTwo.Name, " ")
	//if name1[0] != name2[0] || (string(name1[1][0]) != "1" && string(name2[1][0]) != "2") {
	//	log.Fatalf("Error: Names of fastq pair are not the same\n")
	//}
	curr.Fwd = fqOne
	curr.Rev = fqTwo
	if done1 || done2 {
		return &curr, true
	}
	return &curr, false
}

func ReadFastqsPairs(er *fileio.EasyReader, er2 *fileio.EasyReader) []*PairedEnd {
	var curr *PairedEnd
	var done bool
	var answer []*PairedEnd
	for curr, done = NextFastqPair(er, er2); !done; curr, done = NextFastqPair(er, er2) {
		answer = append(answer, curr)
	}
	return answer
}

func WritePairToFileHandle(file io.Writer, file2 io.Writer, fq []*PairedEnd) {
	for i := 0; i < len(fq); i++ {
		_, err = fmt.Fprintf(file, "%s\n%s\n%s\n%s\n", "@"+fq[i].Fwd.Name, dna.BasesToString(fq[i].Fwd.Seq), "+", string(fq[i].Fwd.Qual))
		common.ExitIfError(err)
		_, err = fmt.Fprintf(file2, "%s\n%s\n%s\n%s\n", "@"+fq[i].Rev.Name, dna.BasesToString(fq[i].Rev.Seq), "+", string(fq[i].Rev.Qual))
		common.ExitIfError(err)
	}
}

func WritePair(readOne string, readTwo string, records []*PairedEnd) {
	file := fileio.EasyCreate(readOne)
	file2 := fileio.EasyCreate(readTwo)
	defer file.Close()
	defer file2.Close()
	WritePairToFileHandle(file, file2, records)
}
