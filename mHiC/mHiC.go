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
	Start1     int
	Start2     int
	Bin1       int
	Bin2       int
	InterChrom bool
	ReadName   string
	Uni        bool
	Full       string
}

type BinContact struct {
	Chrom1 string
	Chrom2 string
	Bin1   int
	Bin2   int
	Count  int
	Key    string
}
type Bias map[string]map[int]float64

type Prior map[int]float64

func ReadPrior(inFile string) Prior {
	var slc []string
	mp := make(map[int]float64)
	in := fileio.Read(inFile)
	for i := range in {
		slc = strings.Split(in[i], "\t")
		mp[parse.StringToInt(slc[0])] = parse.StringToFloat64(slc[1])
	}
	return mp
}

func ReadBias(in string) Bias {
	var val float64
	b := make(map[string]map[int]float64)
	s := fileio.Read(in)
	for i := range s {
		ss := strings.Split(s[i], "\t")
		if parse.StringToFloat64(ss[2]) > 2 || parse.StringToFloat64(ss[2]) < 0.5 {
			val = -1
		} else {
			val = parse.StringToFloat64(ss[2])
		}

		_, found := b[ss[0]]
		if !found {
			b[ss[0]] = make(map[int]float64)
		}
		b[ss[0]][parse.StringToInt(ss[1])] = val
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
	i.Start1 = parse.StringToInt(columns[2])
	i.Chrom2 = columns[6]
	i.Start2 = parse.StringToInt(columns[7])
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

func GoReadBinContactToChan(inFile string) <-chan BinContact {
	var wg sync.WaitGroup
	data := make(chan BinContact, 1000)
	wg.Add(1)
	go readBinContactToChan(inFile, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func readBinContactToChan(inFile string, data chan<- BinContact, wg *sync.WaitGroup) {
	var line BinContact
	var done bool
	in := fileio.EasyOpen(inFile)

	for line, done = nextBinContact(in); !done; line, done = nextBinContact(in) {
		data <- line
	}
	err := in.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func nextBinContact(reader *fileio.EasyReader) (BinContact, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return BinContact{}, true
	}
	return processBinContact(line), false
}

func processBinContact(line string) BinContact {
	var bc BinContact
	slc := strings.Split(line, "\t")
	bc.Chrom1 = slc[0]
	bc.Bin1 = parse.StringToInt(slc[1])
	bc.Chrom2 = slc[2]
	bc.Bin2 = parse.StringToInt(slc[3])
	bc.Count = parse.StringToInt(slc[4])
	bc.Key = strings.Join(slc[0:4], "_")
	return bc
}
