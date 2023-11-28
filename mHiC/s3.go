package mHiC

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
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

type S3settings struct {
	restFragBed      string
	peBam            string
	outDir           string //default is "."
	cutDistanceLow   int    // minimum distance between total read-pair distance and assigned genome cut-site (recommended 50)
	cutDistanceUpper int    // max distance between total read-pair distance and assigned genome cut site (Recommended 500)
	resolution       int    // Size of fixed window or number of RE fragments per bin. Default is 10kb
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

func xaStringToBed(s string) bed.Bed {
	var bd bed.Bed
	fields := strings.Split(s, ",")
	bd.Chrom = fields[0]
	bd.Strand, bd.ChromStart = remove1stChar(fields[1])
	bd.ChromEnd = bd.ChromStart + cigar.ReferenceLength(cigar.FromString(fields[2]))
	bd.FieldsInitialized = 6
	return bd
}

func getXA(a sam.Sam) []interval.Interval {
	var bd bed.Bed
	var outIntervals []interval.Interval
	out, found, _ := sam.QueryTag(a, "XA")
	if !found {
		return []interval.Interval{}
	}
	w := removeTrailingCharAndSplit(out.(string), ";")
	for _, i := range w {
		bd = xaStringToBed(i)
		outIntervals = append(outIntervals, bd)
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

func isSC(a, b interval.Interval) bool {
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

func isRL(frag1, frag2 interval.Interval) bool {
	x, y := orderPairsByPos(frag1, frag2)
	if x.GetChromEnd() == y.GetChromStart() {
		return true
	}
	return false
}

func getInteractionType(aln1, aln2, frag1, frag2 interval.Interval) string {
	var it string
	if frag1 == frag2 { //SC or DE
		if isSC(aln1, aln2) {
			it = "SC"
		} else if isDE(aln1, aln2) {
			it = "DE"
		}
	} else if isRL(frag1, frag2) {
		it = "RL"
	} else {
		it = "VI"
	}
	return it
}

func middleOfInterval(i interval.Interval) interval.Interval {
	mid := i.GetChromStart() + i.GetChromEnd()/2
	bd := bed.Bed{Chrom: i.GetChrom(), ChromStart: mid, ChromEnd: mid + 1, FieldsInitialized: 3}
	return bd
}

func getOverlappingReFrag(tree map[string]*interval.IntervalNode, q interval.Interval) []interval.Interval {
	mid := middleOfInterval(q)
	return interval.Query(tree, mid, "any")
}

func getCisDist(a, b interval.Interval) int {
	if a.GetChrom() != b.GetChrom() {
		return numbers.MaxInt
	}
	return numbers.AbsInt(a.GetChromStart() - b.GetChromStart())
}

func getFragDist(aln1, aln2, frag1, frag2 interval.Interval) int {
	v1 := numbers.Min(numbers.AbsInt(frag1.GetChromEnd()-aln1.GetChromStart()), numbers.AbsInt(aln1.GetChromStart()-frag1.GetChromStart()))
	v2 := numbers.Min(numbers.AbsInt(frag2.GetChromEnd()-aln2.GetChromStart()), numbers.AbsInt(aln2.GetChromStart()-frag2.GetChromStart()))
	return v1 + v2
}

func S3(s S3settings) {
	var a, b sam.Sam
	var pe sam.SamPE
	var alnPos1, alnPos2, sa1, sa2, frag1, frag2 []interval.Interval
	var full bool = true
	var match bool = true
	var i, aln1, aln2, fragDist, cisDist, mid1, mid2 int
	var it string
	var vi []string

	p := getPrefix(s.peBam)                                               // get the filename minus the last extension
	vp := fileio.EasyCreate(fmt.Sprintf("%s/%s.validPairs", s.outDir, p)) //create valid pairs file
	um := fileio.EasyCreate(fmt.Sprintf("%s/%s.umiMulti", s.outDir, p))   //create uni-multi file (I don't know what this file is)
	reTree := createReTree(s.restFragBed)

	inChan, _ := sam.GoReadToChan(s.peBam)

	// sam file MUST be sorted by read-name
	for full {
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
		if a.RName == "" && b.RName == "" { //what is the best way to check for unmapped?
			continue
		}
		pe = sam.SamToPeSamCheckFlag(a, b)
		alnPos1 = append(alnPos1, pe.R1)
		sa1 = getXA(pe.R1)
		for i = range sa1 {
			alnPos1 = append(alnPos1, sa1[i])
		}
		alnPos2 = append(alnPos2, pe.R2)
		sa2 = getXA(pe.R2)
		for i = range sa2 {
			alnPos2 = append(alnPos2, sa2[i])
		}
		for aln1 = range alnPos1 {
			frag1 = getOverlappingReFrag(reTree, alnPos1[aln1])
			if len(frag1) != 1 {
				continue
			}
			for aln2 = range alnPos2 {
				frag2 = getOverlappingReFrag(reTree, alnPos2[aln2])
				if len(frag2) != 1 {
					continue
				}
				it = getInteractionType(alnPos1[aln1], alnPos2[aln2], frag1[0], frag2[0])
				if it == "VI" {
					fragDist = getFragDist(alnPos1[aln1], alnPos2[aln2], frag1[0], frag2[0])
					cisDist = getCisDist(alnPos1[aln1], alnPos2[aln2])
					if fragDist <= s.cutDistanceLow || fragDist >= s.cutDistanceUpper || cisDist <= 2*s.resolution {
						it = "DUMP"
					} else {
						mid1 = 1 //placeholder
						mid2 = 2
						if alnPos1[aln1].GetChrom() == alnPos2[aln2].GetChrom() && mid1 == mid2 {
							it = "DUMP"
						} else {
							vi = append(vi, fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d",
								alnPos1[aln1].GetChrom(), alnPos1[aln1].GetChromStart()))
							//not going to explicitly format by position order, see if it matters later
						}
					}
				}
			}
		}
	}
}
