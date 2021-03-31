package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

// FastqBig is currently used by the aligner.  It holds
// the pre-computed reverse complement of the seq as well
// as the sequences as originally provided.  For both
// of these sequences it also holds
type FastqBig struct {
	Name      string
	Seq       []dna.Base
	SeqRc     []dna.Base
	Qual      []uint8
	Rainbow   []dnaTwoBit.TwoBit
	RainbowRc []dnaTwoBit.TwoBit
}


func ReadBigToChan(filename string, output chan<- FastqBig) {
	var curr Fastq
	var currBig FastqBig
	var done bool
	var err error
	file := fileio.EasyOpen(filename)
	for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
		currBig = ToFastqBig(curr)
		output <- currBig
	}
	close(output)
	err = file.Close()
	exception.PanicOnErr(err)
}

// ToFastqBig converts of Fastq into a FastqBig
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
func ReadFqBig(filename string) (FastqBig, bool) {
	reader := fileio.OpenByteReader(filename)
	answer := FastqBig{}
	line, done := fileio.ReadLine(reader)
	if done {
		return answer, true
	}
	answer.Name = strings.Split(line.String()[1:], " ")[0]
	line, done = fileio.ReadLine(reader)
	if done {
		log.Panicf("Only got a name, but not a sequence at the end of %s\n", filename)
	}
	//set up sequence and reverse comp
	answer.Seq = dna.ByteSliceToDnaBases(line.Bytes())
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)

	//performs two bit conversion
	answer.Rainbow = dnaTwoBit.TwoBitRainbow(answer.Seq)
	answer.RainbowRc = dnaTwoBit.TwoBitRainbow(answer.SeqRc)

	line, done = fileio.ReadLine(reader)
	if done {
		log.Panicf("Error: Only got a name and sequence but no plus sign at the end of %s\n", filename)
	}
	if line.String() != "+" {
		log.Panicf("Error: This line should be a + (plus) sign:\n%s\nIn %s", line.String(), filename)
	}
	line, done = fileio.ReadLine(reader)
	if done {
		log.Panicf("Error: Did not find a set of quality scores at the end of %s\n", filename)
	}
	answer.Qual = ToQual(line.Bytes())
	return answer, false
}
