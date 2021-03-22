package sam

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
	//"github.com/vertgenlab/gonomics/chromInfo"
	"log"
	"os"
	//"sync"
)

//TODO: respond to craigs comments of git
//I am worried about opening all the files, then reading them, then closing them all. What if you put reading the file in the loop where you open it? Then you could: open, read, close, go to next file, open, read close, go to next file, etc
func ReadFilesToChan(filenames []string, output chan<- []SamAln) {
	var curr []SamAln = make([]SamAln, len(filenames))
	var done bool
	var files []*fileio.EasyReader = make([]*fileio.EasyReader, len(filenames))
	for i := 0; i < len(filenames); i++ {
		files[i] = fileio.EasyOpen(filenames[i])
		defer files[i].Close()
		ReadHeader(files[i])
	}
	for curr, done = NextSamSplit(files); done != true; curr, done = NextSamSplit(files) {
		//for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
		output <- curr
	}
	close(output)
}

//TODO: respond to craigs comments
//It is hard to know what this does from the name "ReadFiles", which will be seen as sam.ReadFiles. It seems to read sam files, but also do lots of other stuff
func ReadFiles(filenames []string, output string) {
	os.Remove(output)
	var samfile Sam = Sam{Header: Header{}, Aln: []SamAln{}}
	var curr SamAln
	//var curr []*SamAln = make([]*SamAln, len(filenames))
	var done bool
	var i int
	//nodeMap := make(map[string]uint32)
	var files []*fileio.EasyReader = make([]*fileio.EasyReader, len(filenames))
	samfile.Header = Header{Text: []string{"@HD\tVN:1.6\tSO:Sorted"}, Chroms: nil}
	log.Printf("len=%d", len(filenames))
	for i = 0; i < len(filenames); i++ {
		log.Printf("len=%s", filenames[i])
		files[i] = fileio.EasyOpen(filenames[i])

		next := ReadHeader(files[i])
		samfile.Header.Text = append(samfile.Header.Text, next.Text[1])
		//header.Chroms = append(header.Chroms, next.Chroms[0])

		defer files[i].Close()
	}
	answer := make(map[string]SamAln)
	for _, chr := range files {
		for curr, done = NextAlignment(chr); done != true; curr, done = NextAlignment(chr) {
			answer = getBestSam(curr, answer)
		}
	}
	//nodeId := chromInfo.SliceToMap(samfile.Chroms)
	for key := range answer {
		samfile.Aln = append(samfile.Aln, answer[key])
	}
	SortByCoord(samfile.Aln)
	Write(output, samfile)
	CheckSam(samfile.Aln)
	//
}

func CheckSam(query []SamAln) {
	var yes, no int64 = 0, 0
	for i := 0; i < len(query); i++ {
		if CheckAlignment(query[i]) {
			yes++
			//log.Printf(sam.SamAlnToString(query[i]))
		} else {
			no++
			//log.Printf("This did not map:\n%s\n", sam.SamAlnToString(query[i]))
		}
	}
	log.Printf("Total number of reads aligned: %d...", len(query))
	log.Printf("Number of reads correctly aligned: %d...\n", yes)
	log.Printf("Number of reads mismapped: %d...\n", no)
}

func CheckAlignment(aln SamAln) bool {
	var answer bool = false

	qName := strings.Split(aln.QName, "_")
	refPos := common.StringToInt(qName[1])
	queryPos := getStartRead(aln)
	if queryPos == refPos {
		return true
	}

	return answer
}

func getStartRead(aln SamAln) int {
	var alignedPos int = int(aln.Pos)
	if aln.Cigar[0].Op == 'S' {
		alignedPos += aln.Cigar[0].RunLength
	}
	return alignedPos
}

func getStartNode(aln SamAln) uint32 {
	words := strings.Split(aln.Extra, "\t")
	if !strings.Contains(words[1], "GP:Z:") {
		log.Fatalf("Error: No graph path found...")
	}
	words = strings.Split(words[0][5:], ":")
	return common.StringToUint32(words[0])
}

/*
func FilesChanWorker(input <-chan []*SamAln, output chan<- *SamAln, wg *sync.WaitGroup) {
	//blastzScores := make(map[string]*SamAln)
	for read := range input {
		output <- getBestSam(read)
	}
	wg.Done()
}*/

func NextSamSplit(readers []*fileio.EasyReader) ([]SamAln, bool) {
	var answer []SamAln
	//answer := make([]*SamAln, len(readers))
	var done bool
	var curr SamAln
	for i := 0; i < len(readers)-1; i++ {
		curr, done = NextAlignment(readers[i])
		//answer[i] = curr
		answer = append(answer, curr)
	}
	curr, done = NextAlignment(readers[len(readers)-1])
	answer = append(answer, curr)
	if done == true {
		return answer, true
	}
	return answer, false
}

func getBestSam(curr SamAln, answer map[string]SamAln) map[string]SamAln {
	_, ok := answer[curr.QName]
	if !ok {
		answer[curr.QName] = curr
	} else {
		currBlast := getBlastScore(curr)
		if getBlastScore(answer[curr.QName]) < currBlast {
			answer[curr.QName] = curr
		}
	}
	return answer
}

func getBlastScore(aln SamAln) int64 {
	//var blastz int64 = 0
	words := strings.Split(aln.Extra, "\t")
	if !strings.Contains(words[0], "BZ:i:") {
		log.Fatalf("Error: blastz score was not found...")
	}
	return common.StringToInt64(words[0][5:])

}
