package giraf

import (
	"fmt"
	"github.com/edotau/simpleio"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/fastq"
	"os"
	"sync"
	
	"strconv"
	"bytes"
	"io"
)

type Giraf struct {
	QName     string
	QStart    int
	QEnd      int
	Flag      uint8
	PosStrand bool
	Path      *Path
	Aln       []*cigar.Cigar
	ByteCigar []simpleio.ByteCigar // current cigar will need to be expanded
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

type Path struct {
	TStart int      // The path starts on the TStart base (0-based, closed) of Nodes[0]
	Nodes  []uint32 // The node Id/Index of all the nodes in the path
	TEnd   int      // The path ends on the TEnd base (0-based, open) of Nodes[len(Nodes)-1]
}

type Note struct {
	Tag   string
	Type  rune
	Value string
}

type GirafPair struct {
	Fwd *Giraf
	Rev *Giraf
}

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

func ReadToChan(file *fileio.EasyReader, data chan<- *Giraf, wg *sync.WaitGroup) {
	for curr, done := NextGiraf(file); !done; curr, done = NextGiraf(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

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

func NextGiraf(reader *fileio.EasyReader) (*Giraf, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	return stringToGiraf(line), false
}

func Write(filename string, gfs []*Giraf) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	for i := 0; i < len(gfs); i++ {
		WriteGirafToFileHandle(file, gfs[i])
	}
}

func GirafChanToFile(filename string, input <-chan *Giraf, wg *sync.WaitGroup) {
	file := fileio.MustCreate(filename)
	defer file.Close()

	for line := range input {
		WriteGirafToFileHandle(file, line)
	}
	wg.Done()
}

func GirafPairChanToFile(filename string, input <-chan *GirafPair, wg *sync.WaitGroup) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	for pair := range input {
		WriteGirafToFileHandle(file, pair.Fwd)
		WriteGirafToFileHandle(file, pair.Rev)
	}
	wg.Done()
}

func WriteGirafToFileHandle(file *os.File, gf *Giraf) {
	_, err := fmt.Fprintf(file, "%s\n", GirafToString(gf))
	common.ExitIfError(err)
}

func WriteSimpleGiraf(filename string, input <-chan *Giraf, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer
	
	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	
	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		buf = GirafStringBuilder(gp, buf)
		io.Copy(file, buf)
		buf.Reset()
		simplePool.Put(buf)
	}
	file.Close()
	wg.Done()
}

func GirafStringBuilder(g *Giraf, buf *bytes.Buffer) *bytes.Buffer {
	var err error
	buf = buildGirafString(buf, g, err)
	err = buf.WriteByte('\n')
	common.ExitIfError(err)
	return buf
}

func buildGirafString(buf *bytes.Buffer, g *Giraf, err error) *bytes.Buffer {
	_, err = buf.WriteString(g.QName)
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.QStart))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.QEnd))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(int(g.Flag)))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteRune(common.StrandToRune(g.PosStrand))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(PathToString(g.Path))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(simpleio.ByteCigarString(g.ByteCigar))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.AlnScore))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(int(g.MapQ)))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(simpleio.ByteDnaBasesToString(g.Seq))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(fastq.Uint8QualToString(g.Qual))
	//common.ExitIfError(err)
	//err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(NotesToString(g.Notes))
	common.ExitIfError(err)
	return buf
}
