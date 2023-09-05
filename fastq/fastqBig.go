package fastq

import (
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fileio"
)

// FastqBig holds normal fastq info, as well as the reverse complement of the sequence
// and two rainbow tables where every possible offset of Seq and SeqRc are stored
// in TwoBit format.
type FastqBig struct {
	Name      string
	Seq       []dna.Base
	SeqRc     []dna.Base
	Qual      []uint8
	Rainbow   []dnaTwoBit.TwoBit
	RainbowRc []dnaTwoBit.TwoBit
}

// ReadBigToChan reads a given fastq file into FastqBig structs and sends them into the provided channel.
func ReadBigToChan(filename string, output chan<- FastqBig) {
	var curr Fastq
	var currBig FastqBig
	var done bool
	file := fileio.EasyOpen(filename)
	defer file.Close()
	for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
		currBig = ToFastqBig(curr)
		output <- currBig
	}
	close(output)
}

// ToFastqBig converts a Fastq to FastqBig format.
func ToFastqBig(a Fastq) FastqBig {
	answer := FastqBig{}
	answer.Name = a.Name
	answer.Seq = a.Seq
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)
	answer.Qual = a.Qual
	answer.Rainbow = dnaTwoBit.TwoBitRainbowDeReference(answer.Seq)
	answer.RainbowRc = dnaTwoBit.TwoBitRainbowDeReference(answer.SeqRc)
	return answer
}

// ReadFqBig returns a FastqBig struct that is converted to dnaTwoBit format
// and includes a rainbow offset to perform hash look ups in GSW aligner.
func ReadFqBig(reader *fileio.ByteReader) (FastqBig, bool) {
	answer := FastqBig{}
	line, done := fileio.ReadLine(reader)
	if done {
		return FastqBig{}, true
	}
	answer.Name = strings.Split(line.String()[1:], " ")[0]
	line, done = fileio.ReadLine(reader)
	if done {
		return FastqBig{}, true
	}
	// set up sequence and reverse comp
	answer.Seq = dna.ByteSliceToDnaBases(line.Bytes())
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)

	// performs two bit conversion
	answer.Rainbow = dnaTwoBit.TwoBitRainbowDeReference(answer.Seq)
	answer.RainbowRc = dnaTwoBit.TwoBitRainbowDeReference(answer.SeqRc)

	line, done = fileio.ReadLine(reader)
	if done {
		return FastqBig{}, true
	}
	if line.String() != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	line, done = fileio.ReadLine(reader)
	if done {
		return FastqBig{}, true
	}
	answer.Qual = ToQual(line.Bytes())
	return answer, false
}
