package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
	"sync"
)

// ReadToChan streams the input file to the input data and header channel
// so that only a small portion of the file is kept in memory at a time.
// the header channel will have a single send (the Header struct) then the
// channel will be closed. The header can be retrieved by `header := <-headerChan`.
// Records will continuously stream to the input data channel until the end of the
// file is reached at which point the data channel will be closed.
func ReadToChan(filename string, data chan<- Aln, header chan<- Header) {
	var file *fileio.EasyReader
	var curr Aln
	var done bool
	var err error

	file = fileio.EasyOpen(filename)

	header <- ReadHeader(file)
	close(header)

	for curr, done = ReadNext(file); !done; curr, done = ReadNext(file) {
		data <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(data)
}

// GoReadToChan streams the input file so that only a small portion
// of the file is kept in memory at a time. This function wraps the
// ReadToChan function to automatically handle channel creation,
// goroutine spawning, and header retrieval.
func GoReadToChan(filename string) (<-chan Aln, Header) {
	data := make(chan Aln, 1000)
	header := make(chan Header)
	go ReadToChan(filename, data, header)
	return data, <-header
}

// ReadNext takes a ByteReader and returns the next Sam record as well as a boolean flag
// indicating if the file is finished being read. If there is a Sam record to process
// the function will return the Sam record and 'false'. After processing all Sam records
// in the file, the function will return a blank Aln and 'true'.
func ReadNext(reader *fileio.EasyReader) (Aln, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return Aln{}, true
	}
	return processAlignmentLine(line), false
}

// Read the entire file into a Sam struct where each record
// is an index in Sam.Aln. Note that sam files can get very large
// such that storing the entire file in memory is not feasible. Most
// sam files should be read using the GoReadToChan function which
// streams sam records so only a small portion of the file is kept
// in memory at any given time.
func Read(filename string) ([]Aln, Header) {
	var alignments []Aln
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)

	for curr, doneReading := ReadNext(file); !doneReading; curr, doneReading = ReadNext(file) {
		alignments = append(alignments, curr)
	}

	err := file.Close()
	exception.PanicOnErr(err)
	return alignments, header
}

// processAlignmentLine parses a string representation of a sam file into a single Aln struct.
func processAlignmentLine(line string) Aln {
	var curr Aln
	var err error
	var currUint uint64
	var currInt int64

	words := strings.SplitN(line, "\t", 12)
	if len(words) < 11 {
		log.Fatalf("malformed sam file: was expecting at least 11 columns per line, but this line did not:\n%s\n", line)
	}

	// QName
	curr.QName = words[0]

	// Flag parsing
	currUint, err = strconv.ParseUint(words[1], 10, 16)
	if err != nil {
		log.Fatalf("error processing sam record: could not parse '%s' to a non-negative integer", words[1])
	}
	curr.Flag = uint16(currUint)

	// RName
	curr.RName = words[2]

	// Pos parsing
	currUint, err = strconv.ParseUint(words[3], 10, 32)
	if err != nil {
		log.Fatalf("error processing sam record: could not parse '%s' to a non-negative integer", words[3])
	}
	curr.Pos = uint32(currUint)

	// MapQ parsing
	currUint, err = strconv.ParseUint(words[4], 10, 8)
	if err != nil {
		log.Fatalf("error processing sam record: could not parse '%s' to a non-negative integer", words[4])
	}
	curr.MapQ = uint8(currUint)

	// Cigar parsing
	curr.Cigar = cigar.FromString(words[5])

	// RNext
	curr.RNext = words[6]

	// PNext parsing
	currUint, err = strconv.ParseUint(words[7], 10, 32)
	if err != nil {
		log.Fatalf("error processing sam record: could not parse '%s' to a non-negative integer", words[7])
	}
	curr.PNext = uint32(currUint)

	// TLen parsing
	currInt, err = strconv.ParseInt(words[8], 10, 32)
	if err != nil {
		log.Fatalf("error processing sam record: could not parse '%s' to an integer", words[8])
	}
	curr.TLen = int32(currInt)

	// Seq parsing
	curr.Seq = dna.StringToBases(words[9])

	// Qual parsing
	curr.Qual = words[10]

	// Additional Tag fields
	if len(words) > 11 {
		curr.Extra = words[11]
	}
	return curr
}

// ReadHeaderBytes processes the contiguous header from an EasyReader
// and advances the Reader past the header lines.
func ReadHeader(file *fileio.EasyReader) Header {
	var answer Header
	var done bool
	var line string
	for peek, err := file.Peek(1); err == nil && peek[0] == '@' && !done; peek, err = file.Peek(1) {
		line, done = fileio.EasyNextLine(file)
		answer.Text = append(answer.Text, line)
	}

	answer.Metadata.AllTags, answer.Metadata.Comments = parseTagsAndComments(answer.Text)
	answer.Chroms = getChromInfo(answer.Metadata.AllTags)
	answer.Metadata.Version = getVersion(answer.Metadata.AllTags)
	answer.Metadata.SortOrder = getSortOrder(answer.Metadata.AllTags)
	answer.Metadata.Grouping = getGrouping(answer.Metadata.AllTags)
	return answer
}

// SamChanToFile writes an incoming channel of Aln structs to a file.
func SamChanToFile(incomingSams <-chan Aln, filename string, header Header, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	if header.Text != nil {
		WriteHeaderToFileHandle(file, header)
	}
	for alignedRead := range incomingSams {
		WriteToFileHandle(file, alignedRead)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// TODO: improved header creation + building of additional fields in Header
// FastaHeader converts a set of fasta sequences to a generic Sam Header for write.
func FastaHeader(ref []fasta.Fasta) Header {
	var header Header
	header.Text = append(header.Text, fmt.Sprintf("@HD\t%s\tSO:unsorted", samSpecVersion))
	var words string

	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", ref[i].Name, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, chromInfo.ChromInfo{Name: ref[i].Name, Size: len(ref[i].Seq)})
	}
	return header
}

// TODO: consolidate with FastaHeader?
// ChromInfoSamHeader generates a sam header from input chrominfo.
func ChromInfoSamHeader(chromSize []chromInfo.ChromInfo) Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	for i := 0; i < len(chromSize); i++ {
		header.Text = append(header.Text, makeHeaderRefLine(chromSize[i].Name, chromSize[i].Size))
	}
	return header
}

// makeHeaderRefLine creates an @SQ tagged line to store information about ref sequences
// present in a fasta header.
func makeHeaderRefLine(chromName string, chromSize int) string {
	return fmt.Sprintf("@SQ\tSN:%s\tLN:%d", chromName, chromSize)
}

// Write a header and sam alignments to a file. Note that this requires
// the entire slice of alignments to be present in memory at the same time.
// This can be avoided by repeated calls to WriteToFileHandle on alignments
// retrieved from a channel, such as the output from GoReadToChan.
func Write(filename string, data []Aln, header Header) {
	file := fileio.EasyCreate(filename)
	WriteHeaderToFileHandle(file, header)
	for i := range data {
		WriteToFileHandle(file, data[i])
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes a single Aln struct to the input file.
func WriteToFileHandle(file *fileio.EasyWriter, aln Aln) {
	_, err := fmt.Fprintln(file, ToString(aln))
	exception.PanicOnErr(err)
}

// WriteHeaderToFileHandle writes a sam header to the input file.
func WriteHeaderToFileHandle(file *fileio.EasyWriter, header Header) {
	var err error
	for i := range header.Text {
		_, err = fmt.Fprintln(file, header.Text[i])
		exception.PanicOnErr(err)
	}
}
