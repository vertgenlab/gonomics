package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"path"
	"strings"
)

//different possible interaction types:
// DUMP  -- don't use this read (overlaps 0 or > 1 RE fragments or unknown)
// Singleton Interaction (SI)
// Valid Interaction (VI)
// Self Circle (SC)
// Dangling End (DE)
// Re-ligation (RE)

//output files:
// validPairs --

type pair struct {
	qname   string
	pe      peInterval
	rf1     interval.Interval
	rf2     interval.Interval
	it      string
	bin1    string
	bin2    string
	set     s3settings
	rfTree  map[string]*interval.IntervalNode
	s1      bed.Strand
	s2      bed.Strand
	cisDist int
}

type peInterval struct {
	R1 interval.Interval
	R2 interval.Interval
}

type s3settings struct {
	restFragBed      string
	peBam            string
	outDir           string //default is "."
	cutDistanceLow   int    // minimum distance between total read-pair distance and assigned genome cut-site (recommended read length * 2)
	cutDistanceUpper int    // max distance between total read-pair distance and assigned genome cut site (Recommended 500)
	resolution       int    // Size of fixed window or number of RE fragments per bin. Default is 10kb
	validPairsFile   *fileio.EasyWriter
	uniMultiFile     *fileio.EasyWriter
}

func getPrefix(f string) string {
	file := path.Base(f)
	fields := strings.Split(file, ".")
	prefixSlice := fields[0 : len(fields)-1]
	return strings.Join(prefixSlice, ".")
}

func createReTree(reFile string) map[string]*interval.IntervalNode {
	var inIntervals []interval.Interval
	bd := bed.Read(reFile)
	for i := range bd {
		inIntervals = append(inIntervals, bd[i])
	}
	return interval.BuildTree(inIntervals)
}

func removeTrailingCharAndSplit(s, c string) []string {
	var r []string
	ss := strings.Split(s, c)
	for i := range ss {
		if len(ss[i]) != 0 {
			r = append(r, ss[i])
		}
	}
	return r
}

func remove1stChar(s string) (strand bed.Strand, pos int) {
	for i := range s {
		if i == 0 {
			if s[i] == '+' {
				strand = '+'
			} else if s[i] == '-' {
				strand = '-'
			}
		}
		if i > 0 {
			pos = parse.StringToInt(s[i:])
			break
		}
	}
	return strand, pos
}

func xaStringToSam(s string) sam.Sam {
	var sm sam.Sam
	fields := strings.Split(s, ",")
	sm.RName = fields[0]
	sm.Pos = parse.StringToUint32(fields[1][1:])
	sm.Cigar = cigar.FromString(fields[2])
	if fields[1][0] == '-' {
		sm.Flag = 0x10
	} else if fields[1][0] == '+' {
		sm.Flag = 0
	}
	return sm
}

func getXA(a sam.Sam) []interval.Interval {
	var sm sam.Sam
	var outIntervals []interval.Interval
	out, found, _ := sam.QueryTag(a, "XA")
	if !found {
		return []interval.Interval{}
	}
	w := removeTrailingCharAndSplit(out.(string), ";")
	for _, i := range w {
		sm = xaStringToSam(i)
		outIntervals = append(outIntervals, sm)
	}
	return outIntervals
}

func orderPairsByPos(a, b interval.Interval) (x, y interval.Interval) {
	if a.GetChrom() < b.GetChrom() {
		x = a
		y = b
		return x, y
	} else if a.GetChrom() > b.GetChrom() {
		x = b
		y = a
		return x, y
	}
	if a.GetChromStart() <= b.GetChromStart() {
		x = a
		y = b
		return x, y
	} else {
		x = b
		y = a
		return x, y
	}
}

/*func isSC(a, b interval.Interval) bool {
	x, y := orderPairsByPos(a, b)
	if x.(bed.Bed).Strand == '-' && y.(bed.Bed).Strand == '+' {
		return true
	}
	return false
}

func isDE(a, b interval.Interval) bool {
	x, y := orderPairsByPos(a, b)
	if x.(bed.Bed).Strand == '+' && y.(bed.Bed).Strand == '-' {
		return true
	}
	return false
}
*/

func isRL(frag1, frag2 interval.Interval) bool {
	x, y := orderPairsByPos(frag1, frag2)
	if x.GetChromEnd() == y.GetChromStart() {
		return true
	}
	return false
}

func getInteractionType(pr *pair) {
	pr.it = "DUMP"
	if interval.AreEqual(pr.rf1, pr.rf2) { //SC or DE
		pr.it = "DUMP"
	} else if isRL(pr.rf1, pr.rf2) {
		pr.it = "RL"
	} else {
		pr.it = "VI"
	}
}

func middleOfInterval(i interval.Interval) interval.Interval {
	mid := (i.GetChromStart() + i.GetChromEnd()) / 2
	bd := bed.Bed{Chrom: i.GetChrom(), ChromStart: mid, ChromEnd: mid + 1, FieldsInitialized: 3}
	return bd
}

func getOverlappingReFrag(tree map[string]*interval.IntervalNode, q interval.Interval) []interval.Interval {
	mid := middleOfInterval(q)
	return interval.Query(tree, mid, "any")
}

func getCisDist(pr *pair) {
	if pr.pe.R1.GetChrom() != pr.pe.R2.GetChrom() {
		pr.cisDist = numbers.MaxInt
	} else {
		pr.cisDist = numbers.AbsInt(pr.pe.R1.GetChromStart() - pr.pe.R2.GetChromStart())
	}
}

func getFragDist(pr pair) int {
	v1 := numbers.Min(numbers.AbsInt(pr.rf1.GetChromEnd()-pr.pe.R1.GetChromStart()), numbers.AbsInt(pr.pe.R1.GetChromStart()-pr.rf1.GetChromStart()))
	v2 := numbers.Min(numbers.AbsInt(pr.rf2.GetChromEnd()-pr.pe.R2.GetChromStart()), numbers.AbsInt(pr.pe.R2.GetChromStart()-pr.rf2.GetChromStart()))
	return v1 + v2
}

func multiMapping(p sam.SamPE) bool {
	_, found1, _ := sam.QueryTag(p.R1, "XA")
	_, found2, _ := sam.QueryTag(p.R2, "XA")
	if found1 || found2 {
		return true
	}
	return false
}

func determineInteractionType(pr *pair) {
	r1Frag := getOverlappingReFrag(pr.rfTree, pr.pe.R1)
	r2Frag := getOverlappingReFrag(pr.rfTree, pr.pe.R2)
	if len(r1Frag) != 1 || len(r2Frag) != 1 {
		pr.it = "DUMP"
	} else {
		pr.rf1 = r1Frag[0]
		pr.rf2 = r2Frag[0]
		getInteractionType(pr)
		if pr.it == "VI" {
			filterInteractionType(pr)
		}
	}
}

func filterInteractionType(pr *pair) {
	fragDist := getFragDist(*pr)
	getCisDist(pr)
	if fragDist <= pr.set.cutDistanceLow || fragDist >= pr.set.cutDistanceUpper || pr.cisDist <= 2*pr.set.resolution {
		pr.it = "DUMP"
	} else {
		pr.bin1 = fmt.Sprintf("%s_%d", pr.pe.R1.GetChrom(), pr.pe.R1.GetChromStart()/pr.set.resolution)
		pr.bin2 = fmt.Sprintf("%s_%d", pr.pe.R2.GetChrom(), pr.pe.R2.GetChromStart()/pr.set.resolution)
		if pr.bin1 == pr.bin2 {
			pr.it = "DUMP"
		}
	}
}

func parseMultiMappingFrags(pe sam.SamPE, pr pair) []pair {
	var prSlice []pair
	r1 := []interval.Interval{pe.R1}
	r2 := []interval.Interval{pe.R2}
	_, found1, _ := sam.QueryTag(pe.R1, "XA")
	_, found2, _ := sam.QueryTag(pe.R2, "XA")
	if found1 {
		r1 = append(r1, getXA(pe.R1)...)
	}
	if found2 {
		r2 = append(r2, getXA(pe.R2)...)
	}
	for i := range r1 {
		pr.pe.R1 = r1[i]
		if sam.IsPosStrand(r1[i].(sam.Sam)) {
			pr.s1 = '+'
		} else {
			pr.s1 = '-'
		}
		for j := range r2 {
			pr.pe.R2 = r2[j]
			if sam.IsPosStrand(r2[j].(sam.Sam)) {
				pr.s2 = '+'
			} else {
				pr.s2 = '-'
			}
			determineInteractionType(&pr)
			if pr.it == "VI" {
				prSlice = append(prSlice, pr)
			}
		}
	}
	return prSlice
}

func writeInteraction(pr pair, uniMulti string) {
	var str string
	if pr.cisDist == numbers.MaxInt {
		str = fmt.Sprintf("%s\t%s\t%d\t%c\t%s\t%s\t%s\t%d\t%c\t%s\t%s\t%s\t%s", pr.qname, pr.pe.R1.GetChrom(), pr.pe.R1.GetChromStart(), pr.s1, pr.rf1.(bed.Bed).Name, pr.bin1, pr.pe.R2.GetChrom(), pr.pe.R2.GetChromStart(),
			pr.s2, pr.rf2.(bed.Bed).Name, pr.bin2, "interChrom", uniMulti)
	} else {
		str = fmt.Sprintf("%s\t%s\t%d\t%c\t%s\t%s\t%s\t%d\t%c\t%s\t%s\t%d\t%s", pr.qname, pr.pe.R1.GetChrom(), pr.pe.R1.GetChromStart(), pr.s1, pr.rf1.(bed.Bed).Name, pr.bin1, pr.pe.R2.GetChrom(), pr.pe.R2.GetChromStart(),
			pr.s2, pr.rf2.(bed.Bed).Name, pr.bin2, pr.cisDist, uniMulti)
	}
	fileio.WriteToFileHandle(pr.set.validPairsFile, str)
}

func assignStrandToPair(pe sam.SamPE, pr *pair) {
	if sam.IsPosStrand(pe.R1) {
		pr.s1 = '+'
	} else {
		pr.s1 = '-'
	}
	if sam.IsPosStrand(pe.R2) {
		pr.s2 = '+'
	} else {
		pr.s2 = '-'
	}
}

func characterizePairs(s s3settings) {
	var pr pair
	var validPrs []pair
	var a, b sam.Sam
	var pe sam.SamPE
	var full bool = true
	var match bool = true
	var pass bool
	var c int

	reTree := createReTree(s.restFragBed)

	inChan, _ := sam.GoReadToChan(s.peBam)

	// sam file MUST be sorted by read-name
	for full {
		pr = pair{set: s, rfTree: reTree}
		if !match {
			b = <-inChan
		} else {
			a, full = <-inChan
			b, full = <-inChan
		}
		if a.QName != b.QName {
			b = a
			match = false
			continue
		}
		match = true
		if !full {
			continue
		}
		if a.Cigar[0].Op == '*' && b.Cigar[0].Op == '*' { //what is the best way to check for unmapped?
			continue
		}
		pe, pass = sam.SamToPeSamCheckFlag(a, b)
		if !pass {
			fmt.Println(a)
			fmt.Println(b)
			c++
			break
		}
		pr.qname = a.QName
		if !multiMapping(pe) {
			pr.pe.R1 = pe.R1
			pr.pe.R2 = pe.R2
			assignStrandToPair(pe, &pr)
			determineInteractionType(&pr)
			if pr.it == "VI" {
				writeInteraction(pr, "UNI")
				fileio.WriteToFileHandle(s.uniMultiFile, fmt.Sprintf("%s\t%d", pr.qname, 1))
			}
		} else {
			validPrs = parseMultiMappingFrags(pe, pr)
			switch {
			case len(validPrs) == 1:
				writeInteraction(validPrs[0], "UNI")
				fileio.WriteToFileHandle(s.uniMultiFile, fmt.Sprintf("%s\t%d", pr.qname, 1))
			case len(validPrs) > 1:
				for vp := range validPrs {
					writeInteraction(validPrs[vp], "MULTI")
				}
				fileio.WriteToFileHandle(s.uniMultiFile, fmt.Sprintf("%s\t%d", pr.qname, len(validPrs)))
			default:
				continue
			}
		}
	}
	fmt.Println("reads that couldn't be paired: ", c)
	err := s.uniMultiFile.Close()
	exception.PanicOnErr(err)
	err = s.validPairsFile.Close()
	exception.PanicOnErr(err)
}

func main() {

	flag.Parse()

	p := getPrefix(flag.Arg(0))

	s := s3settings{
		resolution:       10000,
		restFragBed:      flag.Arg(1),
		peBam:            flag.Arg(0),
		outDir:           "testdata/output/",
		cutDistanceLow:   50,
		cutDistanceUpper: 500,
		validPairsFile:   fileio.EasyCreate(fmt.Sprintf("testdata/output/%s.validPairs", p)),
		uniMultiFile:     fileio.EasyCreate(fmt.Sprintf("testdata/output/%s.uniMulti", p)),
	}

	characterizePairs(s)
}
