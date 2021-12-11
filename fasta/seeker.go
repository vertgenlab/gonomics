package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	"os"
)

// Seeker enables random access of fasta sequences using a pre-computed index.
type Seeker struct {
	file readSeekCloser
	idx  Index
}

// Close the Seeker.
func (rs *Seeker) Close() error {
	return rs.file.Close()
}

// TODO
// readSeekCloser is a TEMPORARY github actions is updated as the io.ReadSeekCloser
// was added to the builtin io package in 2020.
type readSeekCloser interface {
	io.Reader
	io.Seeker
	io.Closer
}

// NewSeeker opens a fasta file and an fai index file and enables
// seek functionality so the entire fasta file does not need to be
// present in memory.
//
// If you input an empty string for 'index', NewSeeker tries to find
// the index file as 'fasta'.fai.
func NewSeeker(fasta, index string) *Seeker {
	sr := new(Seeker)
	var err error

	sr.file, err = os.Open(fasta)
	exception.FatalOnErr(err)

	if index == "" {
		index = fasta + ".fai"
	}
	sr.idx = readIndex(index)
	return sr
}

// SeekByName returns a portion of a fasta sequence identified by chromosome name.
// Input start and end should be 0-based start-open end-closed.
func SeekByName(sr *Seeker, chr string, start, end int) []dna.Base {
	idx, ok := sr.idx.nameMap[chr]
	if !ok {
		log.Fatalf("ERROR: could not find sequence for fasta record '%s'\n", chr)
	}
	return seek(sr, sr.idx.chroms[idx], start, end)
}

// SeekByIndex returns a portion of a fasta sequence identified by chromosome index (order in file).
// Input start and end should be 0-based start-open end-closed.
func SeekByIndex(sr *Seeker, chr, start, end int) []dna.Base {
	return seek(sr, sr.idx.chroms[chr], start, end)
}

// seek returns a portion of a fasta sequence retrieved using the input chrOffset.
func seek(sr *Seeker, off chrOffset, start, end int) []dna.Base {
	if start > end || start < 0 {
		log.Panicf("illegal start/end position\n\nstart: %d\nend: %d\n", start, end)
	}
	var startOffset, endOffset int
	startOffset = off.offset + ((start / off.basesPerLine) * off.bytesPerLine) + (start % off.basesPerLine)
	endOffset = off.offset + ((end / off.basesPerLine) * off.bytesPerLine) + (end % off.basesPerLine)

	_, err := sr.file.Seek(int64(startOffset), io.SeekStart)
	exception.PanicOnErr(err)

	data := make([]byte, endOffset-startOffset)
	var bytesRead int
	bytesRead, err = sr.file.Read(data)
	if err != nil && err != io.EOF {
		log.Panic(err)
	}
	data = data[:bytesRead] // in case of read past EOF

	answer := make([]dna.Base, end-start)
	var j int
	for i := range data {
		if data[i] == '>' { // in case of read into next fasta record
			break
		}
		if data[i] == '\r' || data[i] == '\n' {
			continue
		}
		answer[j] = dna.ByteToBase(data[i])
		j++
	}
	return answer[:j] // trim off any extra bases in case EOF
}
