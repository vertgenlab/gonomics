package fasta

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	"os"
)

// Seeker enables random access of fasta sequences using a pre-computed index.
type Seeker struct {
	file readSeekCloser // TODO change to io.ReadSeekCloser after update
	idx  Index
}

// Close the Seeker.
func (rs *Seeker) Close() error {
	return rs.file.Close()
}

// TODO
// readSeekCloser is a TEMPORARY interface until github actions
// is updated to the newer version of Go as the io.ReadSeekCloser
// interface was added to the builtin io package in 2020.
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
func SeekByName(sr *Seeker, chr string, start, end int) ([]dna.Base, error) {
	idx, ok := sr.idx.nameMap[chr]
	if !ok {
		log.Fatalf("ERROR: could not find sequence for fasta record '%s'\n", chr)
	}
	var nextChrStartByte int = -1
	if len(sr.idx.chroms) > idx+1 {
		nextChrStartByte = sr.idx.chroms[idx+1].offset
	}
	return seek(sr, sr.idx.chroms[idx], nextChrStartByte, start, end)
}

// SeekByIndex returns a portion of a fasta sequence identified by chromosome index (order in file).
// Input start and end should be 0-based start-closed end-open.
func SeekByIndex(sr *Seeker, chr, start, end int) ([]dna.Base, error) {
	var nextChrStartByte int = -1
	if len(sr.idx.chroms) > chr+1 {
		nextChrStartByte = sr.idx.chroms[chr+1].offset
	}
	return seek(sr, sr.idx.chroms[chr], nextChrStartByte, start, end)
}

var (
	ErrSeekStartOutsideChr = errors.New("requested start position greater than requested chromosome length, nil output")
	ErrSeekEndOutsideChr   = errors.New("requested bases past end of chr, output truncated")
)

// seek returns a portion of a fasta sequence retrieved using the input chrOffset.
func seek(sr *Seeker, off chrOffset, nextChrStartByte, start, end int) ([]dna.Base, error) {
	var err error
	if start > end || start < 0 {
		log.Panicf("illegal start/end position\n\nstart: %d\nend: %d\n", start, end)
	}
	var startOffset, endOffset int
	startOffset = off.offset + ((start / off.basesPerLine) * off.bytesPerLine) + (start % off.basesPerLine)
	endOffset = off.offset + ((end / off.basesPerLine) * off.bytesPerLine) + (end % off.basesPerLine)

	if nextChrStartByte != -1 && startOffset >= nextChrStartByte { // nextChrStartByte is -1 if request chr is last in file
		return nil, ErrSeekStartOutsideChr
	}

	_, err = sr.file.Seek(int64(startOffset), io.SeekStart)
	exception.PanicOnErr(err)

	data := make([]byte, endOffset-startOffset)
	var bytesRead int
	bytesRead, err = sr.file.Read(data)

	if err != nil { // os.File only returns EOF on subsequent read past EOF
		log.Panic(err)
	}

	if bytesRead < len(data) { // EOF
		data = data[:bytesRead] // in case of read past EOF
		err = ErrSeekEndOutsideChr
	}

	answer := make([]dna.Base, end-start)
	var j int
	for i := range data {
		if data[i] == '>' { // in case of read into next fasta record
			err = ErrSeekEndOutsideChr
			break
		}
		if data[i] == '\r' || data[i] == '\n' {
			continue
		}
		answer[j] = dna.ByteToBase(data[i])
		j++
	}
	return answer[:j], err // trim off any extra bases in case EOF
}
