package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fileio"
)

type FastqBig struct {
	Name      string
	Seq       []dna.Base
	SeqRc     []dna.Base
	Qual      []rune
	Rainbow   []*dnaTwoBit.TwoBit
	RainbowRc []*dnaTwoBit.TwoBit
}

func ReadBigToChan(filename string, output chan<- *FastqBig) {
	var curr *Fastq
	var currBig *FastqBig
	var done bool

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
		currBig = ToFastqBig(curr)
		output <- currBig
	}
	close(output)
}

func ToFastqBig(a *Fastq) *FastqBig {
	answer := &FastqBig{}
	answer.Name = a.Name
	answer.Seq = a.Seq
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)
	answer.Qual = a.Qual
	answer.Rainbow = dnaTwoBit.NewTwoBitRainbow(answer.Seq)
	answer.RainbowRc = dnaTwoBit.NewTwoBitRainbow(answer.SeqRc)
	return answer
}
