package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/mHiC"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"math"
	"sort"
	"strings"
)

type settings struct {
	multiFile  string
	multiKeys  []string
	uniFile    string
	outFile    string
	prior      mHiC.Prior
	resolution int
}

type zstruct struct {
	qnameIndex int
	keyIndex   []int
	alignProb  []float64
	possLoci   int
}

type zrevstruct struct {
	pairIndex int
	qnameList []int
	uniCounts int
}

func bisectLeft(slc []int, val int) int {
	for i := range slc {
		if slc[i] < val {
			continue
		} else {
			return i
		}
	}
	return len(slc)
}

func minInt(a, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}

func getPriorVal(mk string, prior map[int]float64, interProb float64, splineDists []int) float64 {
	slc := strings.Split(mk, "_")
	if slc[0] != slc[2] {
		return interProb
	}
	dist := numbers.AbsInt(parse.StringToInt(slc[3]) - parse.StringToInt(slc[1]))

	idx := minInt(bisectLeft(splineDists, dist), len(splineDists)-1)
	return prior[splineDists[idx]]
}

func readMulti(s settings, interProb float64, splineDists []int) (map[int]map[int]int, map[int][]int, []string, map[string]int, []float64, []float64, map[int]string, int) {
	var pi, prior []float64
	var priorVal float64
	var index int = -1
	var currRead, pk string
	var pairIdx int
	var pairKeys []string

	qnameMap := make(map[int]string)
	z := make(map[int]map[int]int) //qname index -- pairKey index -- 0 (for now)
	zrev := make(map[int][]int)
	keyMap := make(map[string]int)

	for i := range s.multiKeys {
		pairKeys = append(pairKeys, s.multiKeys[i])
		keyMap[s.multiKeys[i]] = i
		priorVal = getPriorVal(s.multiKeys[i], s.prior, interProb, splineDists)
		pi = append(pi, priorVal)
		prior = append(prior, priorVal)
		zrev[i] = []int{}
	}

	// index is qname index
	// pairIdx is key index
	multi := mHiC.GoReadInteractionToChan(s.multiFile)
	for i := range multi {
		pk = strings.Join([]string{i.Chrom1, fileio.IntToString(i.Bin1 * s.resolution), i.Chrom2, fileio.IntToString(i.Bin2 * s.resolution)}, "_")
		if i.ReadName != currRead {
			index++
			qnameMap[index] = i.ReadName
			z[index] = make(map[int]int)
		}
		pairIdx = keyMap[pk]
		z[index][pairIdx] = 0
		zrev[pairIdx] = append(zrev[pairIdx], index)
		currRead = i.ReadName
	}
	return z, zrev, pairKeys, keyMap, pi, prior, qnameMap, index
}

func interChromProb(prior map[int]float64) (float64, []int) {
	var min float64 = math.MaxFloat64
	var splineDists []int
	for i := range prior {
		splineDists = append(splineDists, i)
		if prior[i] > min {
			continue
		}
		min = prior[i]
	}
	sort.Ints(splineDists)
	return min / 2, splineDists
}

func fillZstruct(zstructSlc *[]zstruct, z map[int]map[int]int, qnameMap map[int]string, pairKeys []string) {
	// z is map that is qname index --- pairkey index --- prob
	var pairKeyMap map[int]int    // this makes pairkey -- prob (0)
	var zs zstruct                // one per qname
	for i := 0; i < len(z); i++ { //loop through qname indexes
		zs.qnameIndex = i
		zs.alignProb = []float64{}
		zs.keyIndex = []int{}
		pairKeyMap = z[i]
		zs.possLoci = len(pairKeyMap)
		for j := range pairKeyMap {
			zs.keyIndex = append(zs.keyIndex, j)
			zs.alignProb = append(zs.alignProb, 0)
		}
		*zstructSlc = append(*zstructSlc, zs)
	}
}

func fillZrevstruct(zrevstructSlc *[]zrevstruct, zrev map[int][]int) {
	//zrev is a map that contains pairIndex -- slice of qnameIndexes that want to map there
	var zrs zrevstruct               // one per pairKey
	for i := 0; i < len(zrev); i++ { // loop through pairKeys
		zrs.pairIndex = i
		zrs.qnameList = zrev[i]
		zrs.uniCounts = 0
		*zrevstructSlc = append(*zrevstructSlc, zrs)
	}
}

func readUni(s settings, keyMap map[string]int, zrevstructSlc []zrevstruct) int {
	var found bool
	var index, total int
	bc := mHiC.GoReadBinContactToChan(s.uniFile)
	for i := range bc {
		index, found = keyMap[i.Key]
		if !found {
			continue
		}
		zrevstructSlc[index].uniCounts = i.Count
		total += i.Count
	}
	return total
}

func updateZ(pi, sumPi []float64, zstructSlc []zstruct, threshold float64) int {
	var tmpSum, prob float64
	var pairIdx, selectChange int

	for i := range zstructSlc {
		tmpSum = 0
		for j := range zstructSlc[i].keyIndex {
			tmpSum += pi[zstructSlc[i].keyIndex[j]]
		}
		for j := range zstructSlc[i].keyIndex {
			pairIdx = zstructSlc[i].keyIndex[j]
			prob = pi[pairIdx] / tmpSum
			if (zstructSlc[i].alignProb[j]-threshold)*(prob-threshold) < 0 {
				selectChange++
			}
			zstructSlc[i].alignProb[j] = prob
		}
		sumPi[i] = tmpSum
	}
	return selectChange
}

func updatePi(pi, prior, sumPi []float64, zrevstructSlc []zrevstruct, totalReadN int) float64 {
	var diff, alignProbSum, newPi float64
	for i := range pi {
		alignProbSum = 0
		for j := range zrevstructSlc[i].qnameList {
			alignProbSum += 1.0 / sumPi[zrevstructSlc[i].qnameList[j]]
		}
		newPi = alignProbSum*pi[i] + float64(zrevstructSlc[i].uniCounts) + float64(totalReadN)*prior[i]
		diff += math.Pow(pi[i]-newPi, 2)
		pi[i] = newPi
	}
	return diff
}

func calcProb(pi, prior []float64, zstructSlc []zstruct, zrevstructSlc []zrevstruct, multiReadN, totalReadN int) {
	var threshold, diffThresh, diff float64 = 0.5, 0.01, 0
	var _, selectChangeThresh, selectChange int = 500, 1, 0
	sumPi := make([]float64, multiReadN+1)
	for i := 0; i < 501; i++ {
		fmt.Println("\niteration: ", i)
		selectChange = updateZ(pi, sumPi, zstructSlc, threshold)
		diff = updatePi(pi, prior, sumPi, zrevstructSlc, totalReadN)
		fmt.Printf("selectChange: %d\ndiff: %.12f\n", selectChange, math.Sqrt(diff))
		if math.Sqrt(diff) < diffThresh && selectChange < selectChangeThresh {
			break
		}
	}
}

func writeProb(outFile string, zstructSlc []zstruct, pairKeys []string, qnameMap map[int]string) {
	var key, pos, toWrite string
	o := fileio.EasyCreate(outFile)
	for i := range zstructSlc {
		for j := range zstructSlc[i].keyIndex {
			key = pairKeys[zstructSlc[i].keyIndex[j]]
			pos = strings.Replace(key, "_", "\t", -1)
			toWrite = fmt.Sprintf("%s\t%s\t%.12f", qnameMap[i], pos, zstructSlc[i].alignProb[j])
			fileio.WriteToFileHandle(o, toWrite)
		}
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func assignProb(s settings) {
	var zstructSlc []zstruct
	var zrevstructSlc []zrevstruct

	interProb, splineDists := interChromProb(s.prior)
	z, zrev, pairKeys, keyMap, pi, prior, qnameMap, index := readMulti(s, interProb, splineDists)
	fillZstruct(&zstructSlc, z, qnameMap, pairKeys)
	fillZrevstruct(&zrevstructSlc, zrev)
	totalUni := readUni(s, keyMap, zrevstructSlc)
	calcProb(pi, prior, zstructSlc, zrevstructSlc, index+1, index+totalUni+1)
	writeProb(s.outFile, zstructSlc, pairKeys, qnameMap)
}

func getKeys(multiFile string, resolution int) []string {
	var mk string
	var keys []string
	mp := make(map[string]int)
	m := mHiC.GoReadInteractionToChan(multiFile)
	for i := range m {
		mk = strings.Join([]string{i.Chrom1, fileio.IntToString(i.Bin1 * resolution), i.Chrom2, fileio.IntToString(i.Bin2 * resolution)}, "_")
		mp[mk] = 0
	}
	for i := range mp {
		keys = append(keys, i)
	}
	return keys
}

func main() {
	res := 10000
	var s settings = settings{
		multiFile:  "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/multi.txt",
		multiKeys:  getKeys("/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/multi.txt", res),
		uniFile:    "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/uniContacts.txt",
		outFile:    "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/testOutS6.go.txt",
		prior:      mHiC.ReadPrior("/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/prior.goMHIC"),
		resolution: res,
	}
	fileio.Write("/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/multiKeys.uniq.txt", s.multiKeys)
	assignProb(s)
}
