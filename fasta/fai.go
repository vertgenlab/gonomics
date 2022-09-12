package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

// Index stores the byte offset for each fasta sequencing allowing for efficient random access.
type Index struct {
	chroms  []chrOffset    // for search by index
	nameMap map[string]int // maps chr name to index in chroms
}

// String method for Index enables easy writing with the fmt package.
func (idx Index) String() string {
	answer := new(strings.Builder)
	for i := range idx.chroms {
		answer.WriteString(idx.chroms[i].String())
		answer.WriteByte('\n')
	}
	return answer.String()
}

// chrOffset has offset information about each reference. Equivalent to one line of a fai file.
type chrOffset struct {
	name         string // Name of this reference sequence
	len          int    // Total length of this reference sequence, in bases
	offset       int    // Offset within the FASTA file of this sequence's first base
	basesPerLine int    // The number of bases on each line
	bytesPerLine int    // The number of bytes in each line, including the newline
}

// String method for chrOffset enables easy writing with the fmt package.
func (c chrOffset) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%d", c.name, c.len, c.offset, c.basesPerLine, c.bytesPerLine)
}

// readIndex reads a fai index file to an Index struct that can be used for random access.
func readIndex(filename string) Index {
	file := fileio.EasyOpen(filename)
	var answer Index
	var curr chrOffset
	var line string
	var col []string
	var done bool
	var err error
	for line, done = fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		col = strings.Split(line, "\t")
		if len(col) != 5 {
			log.Fatalf("ERROR: malformed index file: %s\nerror on line:\n%s\n", filename, line)
		}

		curr.name = col[0]
		curr.len, err = strconv.Atoi(col[1])
		exception.PanicOnErr(err)
		curr.offset, err = strconv.Atoi(col[2])
		exception.PanicOnErr(err)
		curr.basesPerLine, err = strconv.Atoi(col[3])
		exception.PanicOnErr(err)
		curr.bytesPerLine, err = strconv.Atoi(col[4])
		exception.PanicOnErr(err)

		answer.chroms = append(answer.chroms, curr)
	}

	err = file.Close()
	exception.PanicOnErr(err)

	answer.nameMap = make(map[string]int)
	for i := range answer.chroms {
		answer.nameMap[answer.chroms[i].name] = i
	}
	return answer
}

// CreateIndex for a fasta file for efficient random access.
func CreateIndex(filename string) Index {
	if strings.HasSuffix(filename, ".gz") {
		log.Fatalf("ERROR: cannot index gzipped file '%s'", filename)
	}
	if !(strings.HasSuffix(filename, ".fa") || strings.HasSuffix(filename, ".fasta")) {
		log.Fatalf("ERROR: '%s' is not a fasta file (.fa or .fasta)", filename)
	}

	file, err := os.Open(filename)
	exception.FatalOnErr(err)
	buf := bufio.NewReader(file)

	var answer Index
	var curr chrOffset
	var line, peek []byte
	var totalBytesRead int

	for err != io.EOF {
		// first, find the name line
		line, err = buf.ReadSlice('\n')
		totalBytesRead += len(line)
		if line[0] != '>' {
			continue
		}
		curr.name = string(bytes.TrimRight(line[1:], "\r\n"))
		curr.offset = totalBytesRead

		// handle corner case of record with empty sequence
		peek, err = buf.Peek(1)
		if err == io.EOF || peek[0] == '>' {
			continue // this excludes the record from the index, which is the same behavior as samtools faidx
		}

		// second, determine the number of bytes per line
		line, err = buf.ReadSlice('\n')
		totalBytesRead += len(line)
		curr.bytesPerLine = len(line)

		// third, determine the number of bases per line
		curr.basesPerLine = len(bytes.TrimRight(line, "\r\n"))

		// fourth, determine the total length of sequence in bases
		curr.len = curr.basesPerLine // for first line we read above
		for peek, err = buf.Peek(1); err != io.EOF && peek[0] != '>'; peek, err = buf.Peek(1) {
			line, err = buf.ReadSlice('\n')
			totalBytesRead += len(line)

			switch {
			case len(line) > curr.bytesPerLine: // more bases than first line
				log.Fatalf("ERROR: cannot index fasta record with different line lengths.\nfile:%s\nrecord:%s\n", filename, curr.name)

			case len(line) < curr.bytesPerLine: // less bases than first line, may be end of seq OR illegal line
				curr.len += curr.basesPerLine - (curr.bytesPerLine - len(line))
				peek, err = buf.Peek(1)
				if err == io.EOF || peek[0] == '>' {
					break
				}
				log.Fatalf("ERROR: cannot index fasta record with different line lengths.\nfile:%s\nrecord:%s\n", filename, curr.name)

			default: // expected line len
				curr.len += curr.basesPerLine
				continue
			}
		}
		answer.chroms = append(answer.chroms, curr)
	}

	answer.nameMap = make(map[string]int)
	for i := range answer.chroms {
		answer.nameMap[answer.chroms[i].name] = i
	}
	return answer
}
