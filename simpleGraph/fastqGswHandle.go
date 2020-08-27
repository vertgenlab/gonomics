package simpleGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"strings"
)

type FastqGsw struct {
	ReadOne fastq.FastqBig
	ReadTwo fastq.FastqBig
}

type GirafGsw struct {
	ReadOne giraf.Giraf
	ReadTwo giraf.Giraf
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- FastqGsw) {
	readOne, readTwo := fileio.NewSimpleReader(fileOne), fileio.NewSimpleReader(fileTwo)
	for fq, done := fqPair(readOne, readTwo); !done; fq, done = fqPair(readOne, readTwo) {
		answer <- *fq
	}
	close(answer)
}

func fqPair(reader1 *fileio.SimpleReader, reader2 *fileio.SimpleReader) (*FastqGsw, bool) {
	var done bool
	var fqOne, fqTwo *fastq.FastqBig = &fastq.FastqBig{}, &fastq.FastqBig{}
	fqOne, done = nextFq(reader1)
	if !done || fqOne != nil {
		fqTwo, done = nextFq(reader2)
		if !done || fqTwo != nil {
			return &FastqGsw{ReadOne: *fqOne, ReadTwo: *fqTwo}, false
		} else {
			log.Fatalf("Error: fastq files do not end at the same time...\n")
		}
	}
	return nil, true
}

func nextFq(reader *fileio.SimpleReader) (*fastq.FastqBig, bool) {
	answer := fastq.FastqBig{}
	line, done := fileio.ReadLine(reader)
	if done {
		return nil, true
	}
	answer.Name = strings.Split(string(line[1:]), " ")[0]
	line, done = fileio.ReadLine(reader)
	if done {
		return nil, true
	}
	//set up sequence and reverse comp
	answer.Seq = ByteSliceToDnaBases(line)
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)

	//performs two bit conversion
	answer.Rainbow = dnaTwoBit.NewTwoBitRainbow(answer.Seq)
	answer.RainbowRc = dnaTwoBit.NewTwoBitRainbow(answer.SeqRc)

	line, done = fileio.ReadLine(reader)
	if done {
		return nil, true
	}
	if string(line) != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}

	line, done = fileio.ReadLine(reader)
	if done {
		return nil, true
	}
	answer.Qual = fastq.ToQualUint8(bytes.Runes(line))
	return &answer, false
}

func NewBigFastqPair() fastq.PairedEndBig {
	return fastq.PairedEndBig{
		Fwd: new(fastq.FastqBig),
		Rev: new(fastq.FastqBig),
	}
}

func ByteToBase(b byte) dna.Base {
	switch b {
	case 'A':
		return dna.A
	case 'C':
		return dna.C
	case 'G':
		return dna.G
	case 'T':
		return dna.T
	case 'N':
		return dna.N
	case 'a':
		return dna.A
	case 'c':
		return dna.C
	case 'g':
		return dna.G
	case 't':
		return dna.T
	case 'n':
		return dna.N
	case '-':
		return dna.Gap
	//VCF uses star to denote a deleted allele
	case '*':
		return dna.Gap
	case '.':
		return dna.Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", b)
		return dna.N
	}
}

func ByteSliceToDnaBases(b []byte) []dna.Base {
	var answer []dna.Base = make([]dna.Base, len(b))
	for i, byteValue := range b {
		answer[i] = ByteToBase(byteValue)
	}
	return answer
}
