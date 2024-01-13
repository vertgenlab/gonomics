package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"sort"
	"strings"
	"sync"
)

type settings struct {
	inFile     string
	chromSizes string
	nBins      int
	resolution int
	outDir     string
}

type contacts struct {
	chrom string
	bin1  int
	bin2  int
	score int
}

type sizeDist struct {
	distance int
	score    int
}

type bin struct {
	index                int
	interactionDistances []int
	nObservedContacts    int
	allDistancesSum      float64
	possiblePairs        int
}

var scaleValue float64 = 1e6

func processContact(line string) contacts {
	var i contacts
	columns := strings.Split(line, "\t")
	i.chrom = columns[0]
	i.bin1 = parse.StringToInt(columns[1])
	i.bin2 = parse.StringToInt(columns[3])
	i.score = parse.StringToInt(columns[4])
	return i
}

func nextContact(reader *fileio.EasyReader) (contacts, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return contacts{}, true
	}
	return processContact(line), false
}

func readContactsToChan(inFile string, data chan<- contacts, wg *sync.WaitGroup) {
	var line contacts
	var done bool
	in := fileio.EasyOpen(inFile)

	for line, done = nextContact(in); !done; line, done = nextContact(in) {
		data <- line
	}
	err := in.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func goReadContactsToChan(inFile string) <-chan contacts {
	var wg sync.WaitGroup
	data := make(chan contacts, 1000)
	wg.Add(1)
	go readContactsToChan(inFile, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func buildMhicPrior(s settings) {
	var total int
	mp := make(map[int]int)
	o := fileio.EasyCreate(fmt.Sprintf("%s/uniSizeDist.txt", s.outDir))
	fileio.WriteToFileHandle(o, "distance\tscore")
	contactChan := goReadContactsToChan(s.inFile)
	for i := range contactChan {
		total += i.score
		mp[numbers.AbsInt(i.bin1-i.bin2)] += i.score
	}
	var dists []sizeDist
	for i := range mp {
		fileio.WriteToFileHandle(o, fmt.Sprintf("%d\t%d", i, mp[i]))
		dists = append(dists, sizeDist{distance: i, score: mp[i]})
	}
	err := o.Close()
	exception.PanicOnErr(err)
	sort.Slice(dists, func(i, j int) bool {
		return dists[i].distance < dists[j].distance
	})
	bins := parseIntoBins(dists, total, s.nBins)
	genomeWidePossible, genomeWideBins := populateBinData(s, bins)
	writeBins(s, bins)
}

func populateBinData(s settings, bins []bin) int {
	var nBin, c, possPairs, currBin, genomeWideTotal, genomeWideBins int
	chrSizes := chromInfo.ReadToSlice(s.chromSizes)
	for i := range chrSizes {
		c = 0
		currBin = 0
		nBin = chrSizes[i].Size/s.resolution + 1
		genomeWideBins += nBin
		for j := 0; j < chrSizes[i].Size+s.resolution; j += s.resolution { // j is interaction distance
			possPairs = nBin - c
			c++
			if j < 2*s.resolution {
				continue
			}
			if j <= bins[currBin].interactionDistances[len(bins[currBin].interactionDistances)-1] {
				bins[currBin].allDistancesSum += (float64(j) / scaleValue) * float64(possPairs)
				bins[currBin].possiblePairs += possPairs
				genomeWideTotal += possPairs
			}
			if j >= bins[currBin].interactionDistances[len(bins[currBin].interactionDistances)-1] {
				currBin++
				if currBin == s.nBins {
					break
				}
			}
		}
	}
	return genomeWideTotal, genomeWideBins
}

func writeBins(s settings, bins []bin) {
	o := fileio.EasyCreate(fmt.Sprintf("%s/binStats.%d.txt", s.outDir, s.nBins))
	fileio.WriteToFileHandle(o, "index\tnumBins\ttotalContacts\tbinStart\tpossiblePairs\tallDistancesSum")
	for i := range bins {
		fileio.WriteToFileHandle(o, fmt.Sprintf("%d\t%d\t%d\t%d\t%d\t%f", i, len(bins[i].interactionDistances), bins[i].nObservedContacts, bins[i].interactionDistances[0], bins[i].possiblePairs, bins[i].allDistancesSum))
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func parseIntoBins(dists []sizeDist, totalInteractions, nBins int) []bin {
	var bins []bin
	var binTotal, currBin, interactionsCounted int
	var distances []int
	interactionsPerBin := totalInteractions / nBins

	for i := range dists {
		interactionsCounted += dists[i].score
		binTotal += dists[i].score
		distances = append(distances, dists[i].distance)
		if binTotal >= interactionsPerBin {
			bins = append(bins, bin{index: currBin, interactionDistances: distances, nObservedContacts: binTotal})
			binTotal = 0
			currBin++
			distances = []int{}
			if currBin < nBins {
				interactionsPerBin = (totalInteractions - interactionsCounted) / (nBins - currBin)
			}
		}
	}
	return bins
}

func main() {
	s := settings{
		inFile:     "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/uniContacts.txt",
		chromSizes: "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/ref/plasmodium.chrom.sizes",
		nBins:      100,
		resolution: 10000,
		outDir:     "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output",
	}
	buildMhicPrior(s)
}
