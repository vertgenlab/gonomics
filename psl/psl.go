// package psl defines the Psl struct containing alignment data as well as functions and methods that operate on the Psl struct
package psl

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

// Psl struct declares data fields and types for Psl struct.
type Psl struct {
	Match       int
	MisMatch    int
	RepeatMatch int
	Ns          int
	QNumIns     int
	QBaseIns    int
	TNumIns     int
	TBaseIns    int
	Strand      string
	QName       string
	QSize       int
	QStart      int
	QEnd        int
	TName       string
	TSize       int
	TStart      int
	TEnd        int
	BlockCount  int
	BlockSize   []int
	QList       []int
	TList       []int
}

// PslReader is a struct that contains a ByteReader with an embedded bufioReader and a 2-D byte slice to reduce memory allocation when reading each line.
type PslReader struct {
	Reader  *fileio.ByteReader
	curr    Psl
	columns []string
	done    bool
}

// NewPslReader will open a given text file and return a pointer to a PslReader struct.
func NewPslReader(filename string) *PslReader {
	return &PslReader{
		Reader: fileio.NewByteReader(filename),
		curr:   Psl{},
	}
}

// Read will process a given text file and parse the data fields of each line and return a slice of Psl structs.
func Read(filename string) []Psl {
	var ans []Psl
	reader := NewPslReader(filename)
	// done is a flag that tells up if we finished reading the psl
	// io.EOF already handled in the fileio.ReadLine() function
	for curr, done := pslLine(reader); !done; curr, done = pslLine(reader) {
		ans = append(ans, *curr)
	}
	return ans
}

// pslLine is a helper function which assigns data fields and return a pointer to a Psl.
func pslLine(reader *PslReader) (*Psl, bool) {
	// done is a flag that tells up if we finished reading the psl
	// io.EOF already handled in the fileio.ReadLine() function
	reader.Reader.Buffer, reader.done = fileio.ReadLine(reader.Reader)
	if !reader.done {
		if reader.Reader.Buffer.Bytes()[0] == '#' {
			// Buffer will advance to the next line if there is a hash tag at the beginning of the line
			reader.Reader.Buffer, reader.done = fileio.ReadLine(reader.Reader)
			if reader.done {
				return nil, true
			}
		}
		reader.columns = strings.Split(reader.Reader.Buffer.String(), "\t")
		if len(reader.columns) == 21 {
			reader.curr.Match = common.StringToInt(reader.columns[0])
			reader.curr.MisMatch = common.StringToInt(reader.columns[1])
			reader.curr.RepeatMatch = common.StringToInt(reader.columns[2])
			reader.curr.Ns = common.StringToInt(reader.columns[3])
			reader.curr.QNumIns = common.StringToInt(reader.columns[4])
			reader.curr.QBaseIns = common.StringToInt(reader.columns[5])
			reader.curr.TNumIns = common.StringToInt(reader.columns[6])
			reader.curr.TBaseIns = common.StringToInt(reader.columns[7])
			reader.curr.Strand = reader.columns[8]
			reader.curr.QName = reader.columns[9]
			reader.curr.QSize = common.StringToInt(reader.columns[10])
			reader.curr.QStart = common.StringToInt(reader.columns[11])
			reader.curr.QEnd = common.StringToInt(reader.columns[12])
			reader.curr.TName = reader.columns[13]
			reader.curr.TSize = common.StringToInt(reader.columns[14])
			reader.curr.TStart = common.StringToInt(reader.columns[15])
			reader.curr.TEnd = common.StringToInt(reader.columns[16])
			reader.curr.BlockCount = common.StringToInt(reader.columns[17])
			reader.curr.BlockSize = fileio.StringToIntSlice(reader.columns[18])
			reader.curr.QList = fileio.StringToIntSlice(reader.columns[19])
			reader.curr.TList = fileio.StringToIntSlice(reader.columns[20])
			return &reader.curr, false
		} else {
			log.Fatalf("Error: number of columns do not equal 21...\n")
			return nil, true
		}
	} else {
		return nil, true
	}
}

// ToString will print a Psl struct to a human readable string.
func ToString(p *Psl) string {
	answer := strings.Builder{}
	answer.WriteString(fileio.IntToString(p.Match))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.MisMatch))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.RepeatMatch))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.Ns))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.QNumIns))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.QBaseIns))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.TNumIns))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.TBaseIns))
	answer.WriteByte('\t')
	answer.WriteString(p.Strand)
	answer.WriteByte('\t')
	answer.WriteString(p.QName)
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.QSize))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.QStart))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.QEnd))
	answer.WriteByte('\t')
	answer.WriteString(p.TName)
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.TSize))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.TStart))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.TEnd))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntToString(p.BlockCount))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntSliceToString(p.BlockSize))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntSliceToString(p.QList))
	answer.WriteByte('\t')
	answer.WriteString(fileio.IntSliceToString(p.TList))
	return answer.String()
}
