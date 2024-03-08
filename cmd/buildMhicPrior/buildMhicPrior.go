package main

import (
	"fmt"
	"github.com/cnkei/gospline"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/mHiC"
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
	bias       mHiC.Bias
}

type contacts struct {
	chrom string
	bin1  int
	bin2  int
	score int
	dist  int
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
	i.dist = numbers.AbsInt(parse.StringToInt(columns[3]) - parse.StringToInt(columns[1]))
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
	populateBinData(s, bins)
	//genomeWideBins, maxPossDist, possIntraInRange, possInterAllCount, interChromProb, baselineIntraChromProb := populateBinData(s, bins)
	x, y := calculateProbabilities(s, bins, total) // x = avg distance, y = avg contacts
	writeBins(s, bins)
	splineFit(x, y, s, mp, total)
}

func populateBinData(s settings, bins []bin) (int, int, int, int, float64, float64) {
	var nBin, c, possPairs, currBin, genomeWideIntraTotal, genomeWideBins, possInterAllCount, possIntraInRange, maxPossDist int
	var possIntraAllCount float64
	chrSizes := chromInfo.ReadToSlice(s.chromSizes)

	for i := range chrSizes {
		c = 0
		currBin = 0
		nBin = chrSizes[i].Size/s.resolution + 1
		genomeWideBins += nBin
		for j := 0; j < chrSizes[i].Size+s.resolution; j += s.resolution { // j is interaction distance
			maxPossDist = numbers.Max(maxPossDist, j)
			possPairs = nBin - c
			c++
			if j < s.resolution {
				continue
			}
			possIntraInRange += possPairs
			//fmt.Println(possIntraInRange)
			if j <= bins[currBin].interactionDistances[len(bins[currBin].interactionDistances)-1] {
				bins[currBin].allDistancesSum += (float64(j) / scaleValue) * float64(possPairs)
				bins[currBin].possiblePairs += possPairs
				genomeWideIntraTotal += possPairs
			}
			if j >= bins[currBin].interactionDistances[len(bins[currBin].interactionDistances)-1] {
				currBin++
				if currBin == s.nBins {
					break
				}
			}
		}
		possInterAllCount += nBin * (genomeWideBins - nBin)
		possIntraAllCount += float64(nBin*(nBin+1)) / 2
	}
	interChromProb := 1 / float64(possInterAllCount)
	baselineIntraChromProb := 1 / possIntraAllCount

	return genomeWideBins, maxPossDist, possIntraInRange, possInterAllCount, interChromProb, baselineIntraChromProb
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

func calculateProbabilities(s settings, bins []bin, total int) (avgDistance, avgContacts []float64) {
	var avgC, avgD float64
	o := fileio.EasyCreate(fmt.Sprintf("%s/averageContactsByDist.txt", s.outDir))
	fileio.WriteToFileHandle(o, "averageDistance\taverageContacts")
	for i := range bins {
		avgC = float64(float64(bins[i].nObservedContacts)/float64(bins[i].possiblePairs)) / float64(total) //possible pairs is a bit of normalizing factor
		avgD = (bins[i].allDistancesSum / float64(bins[i].possiblePairs)) * 1000000                        // take closer look at this
		avgContacts = append(avgContacts, avgC)
		avgDistance = append(avgDistance, avgD)
		fileio.WriteToFileHandle(o, fmt.Sprintf("%f\t%.15f", avgD, avgC))
	}
	err := o.Close()
	exception.PanicOnErr(err)
	return avgDistance, avgContacts
}

func splineFit(x, y []float64, s settings, mp map[int]int, total int) {
	var vals, spVals []float64
	o := fileio.EasyCreate(fmt.Sprintf("%s/splineVals.txt", s.outDir))
	fileio.WriteToFileHandle(o, "dist\tsplineVal")
	for i := range mp {
		if i >= 2*s.resolution && float64(i) <= x[len(x)-1] {
			vals = append(vals, float64(i))
		}
	}
	sort.Float64s(vals)
	sp := gospline.NewCubicSpline(x, y)
	for i := range vals {
		spVals = append(spVals, sp.At(vals[i]))
	}
	pValues(s, vals, spVals, x, total)
}

func pValues(s settings, splineX []float64, splineY []float64, avgDist []float64, total int) {
	var b1, b2, dist, priorP, p, place float64
	var pVals [][]float64
	var idx int
	var min, max float64 = minFloatSlice(avgDist), maxFloatSlice(avgDist)
	ch := goReadContactsToChan(s.inFile)
	for i := range ch {
		b1 = s.bias[i.chrom][i.bin1]
		b2 = s.bias[i.chrom][i.bin2]
		if b1 <= 0 || b2 <= 0 {
			place++
			continue
		}
		if float64(i.dist) < min {
			dist = min
		} else if float64(i.dist) > max {
			dist = max
		} else {
			dist = float64(i.dist)
		}
		idx = numbers.Min(bisectLeft(splineX, dist), len(avgDist)-1)
		priorP = splineY[idx] * b1 * b2
		p = numbers.BinomialLeftSummation(total, i.score-1, priorP, false)
		if i.chrom == "chr13" && i.bin1 == 1250000 && i.bin2 == 1230000 {
			fmt.Println(total, i.score-1, priorP, 1-p)
		}
		pVals = append(pVals, []float64{place, 1 - p})
		place++
	}
	sortPValsRank(pVals)
	pValsBH := numbers.BenjaminiHochberg(pVals)
	sortPValsIndex(pValsBH)
	sortPValsIndex(pVals)
	o := fileio.EasyCreate(s.outDir + "uniContactsWithPVals.txt")
	ch = goReadContactsToChan(s.inFile)
	var c int
	idx = 0
	for i := range ch {
		if c > len(pVals)-1 {
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d\t%d\t%d\tnan\tnan", i.chrom, i.bin1, i.bin2, i.score))
			continue
		}
		if float64(idx) == pVals[c][0] {
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d\t%d\t%d\t%.12f\t%.12f", i.chrom, i.bin1, i.bin2, i.score, pVals[c][1], pValsBH[c][1]))
			c++
		} else {
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d\t%d\t%d\tnan\tnan", i.chrom, i.bin1, i.bin2, i.score))
		}
		idx++
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func sortPValsRank(pVals [][]float64) {
	sort.Slice(pVals, func(i, j int) bool {
		return pVals[i][1] < pVals[j][1]
	})
}

func sortPValsIndex(pVals [][]float64) {
	sort.Slice(pVals, func(i, j int) bool {
		return pVals[i][0] < pVals[j][0]
	})
}

func bisectLeft(slc []float64, val float64) int {
	for i := range slc {
		if slc[i] < val {
			continue
		} else {
			return i
		}
	}
	return len(slc)
}

func minFloatSlice(s []float64) float64 {
	f := s[0]
	for i := range s {
		if s[i] < f {
			f = s[i]
		}
	}
	return f
}

func maxFloatSlice(s []float64) float64 {
	f := s[0]
	for i := range s {
		if s[i] > f {
			f = s[i]
		}
	}
	return f
}

func main() {
	s := settings{
		inFile:     "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/uniContacts.txt",
		chromSizes: "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/ref/plasmodium.chrom.sizes",
		nBins:      100,
		resolution: 10000,
		outDir:     "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/",
		bias:       mHiC.ReadBias("/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/biasFile.txt"),
	}
	buildMhicPrior(s)
}
