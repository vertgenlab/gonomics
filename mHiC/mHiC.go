package mHiC

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
	"sync"
)

type Interaction struct {
	Chrom1     string
	Chrom2     string
	Bin1       int
	Bin2       int
	InterChrom bool
	ReadName   string
	Uni        bool
	Full       string
}

type Bias struct {
	chrom string
	bin   int
	value float64
}

func ReadBias(in string) []Bias {
	var b []Bias
	var val float64
	s := fileio.Read(in)
	for i := range s {
		ss := strings.Split(s[i], "\t")
		if parse.StringToFloat64(ss[2]) > 2 || parse.StringToFloat64(ss[2]) < 0.5 {
			val = -1
		} else {
			val = parse.StringToFloat64(ss[2])
		}
		b = append(b, Bias{
			chrom: ss[0],
			bin:   parse.StringToInt(ss[1]),
			value: val,
		})
	}
	return b
}

func readInteractionToChan(inFile string, data chan<- Interaction, wg *sync.WaitGroup) {
	var line Interaction
	var done bool
	in := fileio.EasyOpen(inFile)

	for line, done = nextInteraction(in); !done; line, done = nextInteraction(in) {
		data <- line
	}
	err := in.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func GoReadInteractionToChan(inFile string) <-chan Interaction {
	var wg sync.WaitGroup
	data := make(chan Interaction, 1000)
	wg.Add(1)
	go readInteractionToChan(inFile, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func nextInteraction(reader *fileio.EasyReader) (Interaction, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Interaction{}, true
	}
	return processInteraction(line), false
}

func processInteraction(line string) Interaction {
	var i Interaction
	i.Full = line
	columns := strings.Split(line, "\t")
	i.ReadName = columns[0]
	i.Chrom1 = columns[1]
	i.Chrom2 = columns[6]
	i.Bin1 = binToInt(columns[5])
	i.Bin2 = binToInt(columns[10])
	if columns[12] == "UNI" {
		i.Uni = true
	} else {
		i.Uni = false
	}
	if columns[11] == "interChrom" {
		i.InterChrom = true
	} else {
		i.InterChrom = false
	}
	return i
}

func binToInt(bin string) int {
	columns := strings.Split(bin, "_")
	return parse.StringToInt(columns[1])
}
