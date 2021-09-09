package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strconv"
	"strings"
	"sync"
)

// readSamToChan streams the input Sam file to the input data and header channel
// so that only a small portion of the file is kept in memory at a time.
// the header channel will have a single send (the Header struct) then the
// channel will be closed. The header can be retrieved by `header := <-headerChan`.
// Records will continuously stream to the input data channel until the end of the
// file is reached at which point the data channel will be closed.
func readSamToChan(filename string, data chan<- Sam, header chan<- Header) {
	var file *fileio.EasyReader
	var curr Sam
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

// readSamToChanRecycle is similar to readSamToChan however it reuses Sam structs from receiveRecords avoiding
// allocations for each new read as is necessary in readSamToChan.
func readSamToChanRecycle(filename string, sendRecords chan<- *Sam, receiveRecords <-chan *Sam, header chan<- Header) {
	var file *fileio.EasyReader
	var curr *Sam
	var done bool
	var err error

	file = fileio.EasyOpen(filename)

	header <- ReadHeader(file)
	close(header)

	for curr = range receiveRecords {
		done = readNextRecycle(file, curr)
		if done {
			break
		}
		sendRecords <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(sendRecords)
}

// readBamToChan streams the input Bam file to the input data and header channel
// so that only a small portion of the file is kept in memory at a time.
// the header channel will have a single send (the Header struct) then the
// channel will be closed. The header can be retrieved by `header := <-headerChan`.
// Records will continuously stream to the input data channel until the end of the
// file is reached at which point the data channel will be closed.
func readBamToChan(filename string, data chan<- Sam, header chan<- Header) {
	var file *BamReader
	var head Header
	var err error

	file, head = OpenBam(filename)

	header <- head
	close(header)

	for {
		var curr Sam
		_, err = DecodeBam(file, &curr)
		if err == io.EOF {
			break
		}
		data <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(data)
}

// readBamToChanRecycle is similar to readBamToChan however it reuses Sam structs from receiveRecords avoiding
// allocations for each new read as is necessary in readBamToChan.
func readBamToChanRecycle(filename string, sendRecords chan<- *Sam, receiveRecords <-chan *Sam, header chan<- Header) {
	var file *BamReader
	var head Header
	var curr *Sam
	var err error

	file, head = OpenBam(filename)

	header <- head
	close(header)

	for curr = range receiveRecords {
		_, err = DecodeBam(file, curr)
		if err == io.EOF {
			break
		}
		sendRecords <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(sendRecords)
}

// GoReadToChan streams the input file so that only a small portion
// of the file is kept in memory at a time. This function automatically
// handles channel creation, goroutine spawning, and header retrieval.
//
// GoReadToChan will detect if the input file ends in ".bam" and
// automatically switch to a bam parser.
func GoReadToChan(filename string) (<-chan Sam, Header) {
	data := make(chan Sam, 1000)
	header := make(chan Header)

	if strings.HasSuffix(filename, ".bam") {
		go readBamToChan(filename, data, header)
	} else {
		go readSamToChan(filename, data, header)
	}
	return data, <-header
}

// GoReadToChanRecycle streams the input file so that only a small portion
// of the file is kept in memory at a time. This function automatically
// handles channel creation, goroutine spawning, and header retrieval.
//
// GoReadToChanRecycle will detect if the input file ends in ".bam" and
// automatically switch to a bam parser.
//
// Unlike GoReadToChan, GoReadToChanRecycle has 2 channel returns, a
// receiver channel (parsedRecords) and a sender channel (recycledStructs).
// Parsed sam records can be accessed from the receiver channel. Once the
// received struct is no longer needed, it can be returned to this function
// via the sender channel. Compared to GoReadToChan, this system significantly
// reduces memory allocations as structs can be reused rather then discarded.
// The total number of Sam structs that are allocated can be set with the
// bufferSize input variable.
//
// Note that if bufferSize number of structs are retained by the calling function
// and not returned to this function. No new records will be sent and the
// goroutines will deadlock.
//
// One notable disadvantage is that each read is only valid until the struct
// is returned to this function. This means that any information that must be
// retained between reads must be copied by the calling function.
func GoReadToChanRecycle(filename string, bufferSize int) (parsedRecords <-chan *Sam, recycledStructs chan<- *Sam, header Header) {
	parsedRecordsInit := make(chan *Sam, bufferSize)
	recycledStructsInit := make(chan *Sam, bufferSize)
	headerChan := make(chan Header)

	if strings.HasSuffix(filename, ".bam") {
		go readBamToChanRecycle(filename, parsedRecordsInit, recycledStructsInit, headerChan)
	} else {
		go readSamToChanRecycle(filename, parsedRecordsInit, recycledStructsInit, headerChan)
	}

	// allocate Sam structs and send them to the buffer
	for i := 0; i < bufferSize; i++ {
		recycledStructsInit <- new(Sam)
	}

	// these values are initiated as a seperate variable so that we can have
	// both named return values and send/receive protected channels as the
	// channels need to be input to the readToChan functions with reverse
	// polarity as the returned channels.
	parsedRecords = parsedRecordsInit
	recycledStructs = recycledStructsInit
	header = <-headerChan
	return
}

// ReadNext takes an EasyReader and returns the next Sam record as well as a boolean flag
// indicating if the file is finished being read. If there is a Sam record to process
// the function will return the Sam record and 'false'. After processing all Sam records
// in the file, the function will return a blank Sam and 'true'. ReadNext will advance
// the reader past all header lines beginning with '@'.
func ReadNext(reader *fileio.EasyReader) (Sam, bool) {
	var answer Sam
	var line string
	var done bool

	// read first non-header line
	for line, done = fileio.EasyNextLine(reader); !done && line[0] == '@'; line, done = fileio.EasyNextLine(reader) {
	}

	if done {
		return Sam{}, true
	}

	answer = processAlignmentLine(line)
	return answer, done
}

// readNextRecycle functions similarly to ReadNext, but reuses a Sam struct to reduce
// memory allocations.
func readNextRecycle(reader *fileio.EasyReader, s *Sam) bool {
	var line string
	var done bool

	// read first non-header line
	for line, done = fileio.EasyNextLine(reader); !done && line[0] == '@'; line, done = fileio.EasyNextLine(reader) {
	}

	if done {
		return true
	}

	processAlignmentLineRecycle(line, s)
	return done
}

// Read the entire file into a Sam struct where each record
// is an index in Sam.Sam. Note that sam files can get very large
// such that storing the entire file in memory is not feasible. Most
// sam files should be read using the GoReadToChan function which
// streams sam records so only a small portion of the file is kept
// in memory at any given time.
func Read(filename string) ([]Sam, Header) {
	var alignments []Sam
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)

	for curr, doneReading := ReadNext(file); !doneReading; curr, doneReading = ReadNext(file) {
		alignments = append(alignments, curr)
	}

	err := file.Close()
	exception.PanicOnErr(err)
	return alignments, header
}

// processAlignmentLine parses a string representation of a sam file into a single Sam struct.
func processAlignmentLine(line string) Sam {
	var answer Sam
	processAlignmentLineRecycle(line, &answer)
	return answer
}

// processAlignmentLineRecycle parses a string representation of a sam file into a single Sam struct.
func processAlignmentLineRecycle(line string, curr *Sam) {
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
	curr.Seq = dna.StringToBases(words[9]) // TODO mem efficient dna parsing

	// Qual parsing
	curr.Qual = words[10]

	// Additional Tag fields
	if len(words) > 11 {
		curr.Extra = words[11]
	} else {
		curr.Extra = ""
	}
}

// ReadHeader processes the contiguous header from an EasyReader
// and advances the Reader past the header lines.
func ReadHeader(file *fileio.EasyReader) Header {
	var answer Header
	var done bool
	var line string
	var err error
	var peek []byte

	for peek, err = file.Peek(1); err == nil && peek[0] == '@' && !done; peek, err = file.Peek(1) {
		line, done = fileio.EasyNextLine(file)
		answer.Text = append(answer.Text, line)
	}

	if err != nil && err != io.EOF {
		log.Fatalf("Error: had the following problem while reading the sam header: %s\n", err)
	}

	answer = ParseHeaderText(answer)
	return answer
}

// WriteFromChan writes an incoming channel of Sam structs to a file.
// The input wait group will be decremented once the write finishes.
// This is necessary to ensure the main thread will not terminate
// before the write has finished. Note that the function sending Sam
// to the data channel and WriteFromChan must not be run on the same
// thread, or else the process will deadlock.
func WriteFromChan(data <-chan Sam, filename string, header Header, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	if header.Text != nil {
		WriteHeaderToFileHandle(file, header)
	}
	for record := range data {
		WriteToFileHandle(file, record)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GenerateHeader creates a header given some information about the input sam.
// The chromSize can be either from a chromsizes file or by calling fasta.ToChromInfo
// on a []Fasta. The sortOrder and groupings are defined in metadata.go. Often
// the most relevant are sam.Unsorted (for sortOrder) and sam.None (for grouping) respectively.
// Additional tag information should be generated prior to calling this function and
// passed as a []string which will be appended to the raw text, and parsed into a HeaderTagMap.
func GenerateHeader(chromSize []chromInfo.ChromInfo, additional []string, sortOrder SortOrder, grouping Grouping) Header {
	var header Header
	header.Text = append(header.Text, fmt.Sprintf("@HD\tVN:%s\tSO:%s", samSpecVersion, sortOrder))
	if grouping != None {
		header.Text[0] += fmt.Sprintf("\tGO:%s", grouping)
	}
	for i := 0; i < len(chromSize); i++ {
		header.Text = append(header.Text, makeHeaderRefLine(chromSize[i].Name, chromSize[i].Size))
	}
	header.Text = append(header.Text, additional...)
	return ParseHeaderText(header)
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
func Write(filename string, data []Sam, header Header) {
	file := fileio.EasyCreate(filename)
	if strings.HasSuffix(filename, ".bam") {
		wr := NewBamWriter(file, header)
		for i := range data {
			WriteToBamFileHandle(wr, data[i], 0)
		}
	} else {
		WriteHeaderToFileHandle(file, header)
		for i := range data {
			WriteToFileHandle(file, data[i])
		}
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes a single Sam struct to the input file.
func WriteToFileHandle(file io.Writer, aln Sam) {
	if bamWriter, ok := file.(*BamWriter); ok {
		WriteToBamFileHandle(bamWriter, aln, 0)
	} else {
		_, err := fmt.Fprintln(file, ToString(aln))
		exception.PanicOnErr(err)
	}
}

// WriteHeaderToFileHandle writes a sam header to the input file.
func WriteHeaderToFileHandle(file io.Writer, header Header) {
	if _, ok := file.(*BamWriter); ok {
		return//for bam files, the header is written by WriteToFileHandle already, so we can ignore manual calls to the writeHeaderToFileHandle function.
		//this is a bit messy, excited to hear other people's thoughts at code review
	}
	var err error
	for i := range header.Text {
		_, err = fmt.Fprintln(file, header.Text[i])
		exception.PanicOnErr(err)
	}
}
