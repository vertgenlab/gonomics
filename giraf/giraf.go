// giraf package contains utilities and software to operate genome graph alignments in giraf format
package giraf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"sync"
)

// Giraf struct contains data fields that describes query sequences aligned to a genome graph reference
type Giraf struct {
	QName     string
	QStart    int
	QEnd      int
	Flag      uint8
	PosStrand bool
	Path      Path
	Cigar     []cigar.ByteCigar
	AlnScore  int
	MapQ      uint8
	Seq       []dna.Base // dnaTwoBit?
	Qual      []uint8
	Notes     []Note // Similar to sam, this is should be a list of notes.
	// Each note should be of the form TAG:TYPE:VALUE
	// TAG is two characters
	// TYPE is a single character
	// VALUE will be stored as a string and can then be de-coded based on type
	// An example would be "BZ:i:4000
}

// Path describes an alignment using a traversal slice of Node Ids as well as a start and end coordinate
type Path struct {
	TStart int      // The path starts on the TStart base (0-based, closed) of Nodes[0]
	Nodes  []uint32 // The node Id/Index of all the nodes in the path
	TEnd   int      // The path ends on the TEnd base (0-based, open) of Nodes[len(Nodes)-1]
}

// Note is a scruct containing additional information from the graph alignment
type Note struct {
	Tag   []byte
	Type  byte
	Value string
}

// GirafPair is a paired end giraf alignment with pointers to
// forward (readOne alignment) and reverse (readTwo alignment)
type GirafPair struct {
	Fwd Giraf
	Rev Giraf
}

// Read will process a text file, parse the data fields and assign values to the appropriate giraf data fields
func Read(filename string) []*Giraf {
	var answer []*Giraf
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var curr *Giraf
	var done bool
	for curr, done = NextGiraf(file); !done; curr, done = NextGiraf(file) {
		answer = append(answer, curr)
	}
	return answer
}

// ReadToChan will read a giraf format file to a goroutine channel containing giraf pointers
func ReadToChan(file *fileio.EasyReader, data chan<- *Giraf, wg *sync.WaitGroup) {
	for curr, done := NextGiraf(file); !done; curr, done = NextGiraf(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

// GoReadToChan is a wrapper function which will set up a wait group to ensure the giraf channel is close properly.
// User does not need to worry about closing the channel when process is finished
func GoReadToChan(filename string) <-chan *Giraf {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan *Giraf)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

// NextGiraf performs (*bufio.Reader)ReadString("\n") to buffer reading one line at a time.
// In addition, the function will perform error handling and return a bool == true when EOF is reached.
func NextGiraf(reader *fileio.EasyReader) (*Giraf, bool) {
	line, done := fileio.EasyNextLine(reader)
	if !done {
		return stringToGiraf(line), false
	} else {
		return nil, true
	}
}

// Write will write a slice of giraf structs to a file and perform error handling.
func Write(filename string, gfs []*Giraf) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	for i := 0; i < len(gfs); i++ {
		WriteGirafToFileHandle(file, gfs[i])
	}
}

// GirafChanToFile will write a giraf channel to a file
func GirafChanToFile(filename string, input <-chan *Giraf, wg *sync.WaitGroup) {
	file := fileio.MustCreate(filename)
	defer file.Close()

	for line := range input {
		WriteGirafToFileHandle(file, line)
	}
	wg.Done()
}

// GirafPairChanToFile will write a giraf pair end alignment channel to a file
func GirafPairChanToFile(filename string, input <-chan GirafPair, wg *sync.WaitGroup) {
	var err error
	file := fileio.EasyCreate(filename)
	for pair := range input {
		WriteGirafToFileHandle(file, &pair.Fwd)
		WriteGirafToFileHandle(file, &pair.Rev)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// WriteGirafToFileHandle will write giraf to a io.Writer
func WriteGirafToFileHandle(file io.Writer, gf *Giraf) {
	_, err := fmt.Fprintf(file, "%s\n", GirafToString(gf))
	common.ExitIfError(err)
}
