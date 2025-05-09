package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

type BulkOutputSeqSettings struct {
	InSam       string
	InBed       string
	OutCounts   string
	SingleBx    bool
	DualBx      bool
	InputNorm   string
	PairedEnd   bool
	RepNum      int
	ZScore      string
	CheckBx     string
	IntronStats bool
	BxHopping   string
}

type construct struct {
	name          string
	rawCounts     float64
	normCounts    float64
	normFactor    float64
	intronCount   float64
	intronPercent float64
	zscore        float64
	barcode       []dna.Base
	family        string
	halfCounts    int
}

func processPairEndSam(s BulkOutputSeqSettings, tree map[string]*interval.IntervalNode, countsMap map[string]construct, ref []fasta.Fasta) []string {
	var ansR1, ansR2 []interval.Interval
	var totPair, qualIntronFilter, twoCollapseToOne, same, bothmultianddiff, none, other, totLen3_save, totLen3_discard, pairsFiltered, oneOnlyOtherFilterd, toAdd, bxMismatch, twoBxHalfCounts, oneBxHalfCount, oneOnlyOtherUnmap int
	var nameR1, nameR2 string
	var intron bool
	var curr construct

	peReads, _ := sam.GoReadSamPeToChan(s.InSam)

	if s.BxHopping != "" {
		readBxHoppingFile(countsMap, s.BxHopping)
	}

	for pair := range peReads {
		totPair++

		if pair.R1.MapQ == 0 {
			qualIntronFilter++
			pair.R1.RName = ""
		} else if !filterIntrons(pair.R1) {
			qualIntronFilter++
			pair.R1.RName = ""
		}
		if pair.R2.MapQ == 0 {
			qualIntronFilter++
			pair.R2.RName = ""
		} else if !filterIntrons(pair.R2) {
			qualIntronFilter++
			pair.R2.RName = ""
		}
		if pair.R1.RName == "" && pair.R2.RName == "" {
			pairsFiltered++
			continue
		}

		intron = false
		if detectIntron(pair.R1) || detectIntron(pair.R2) {
			intron = true
		}

		if s.DualBx || s.SingleBx {
			ansR1, toAdd = bxOverlap(pair.R1, tree, ref, s, countsMap)
			bxMismatch += toAdd
			ansR2, toAdd = bxOverlap(pair.R2, tree, ref, s, countsMap)
			bxMismatch += toAdd
		} else {
			ansR1 = interval.Query(tree, pair.R1, "any")
			ansR2 = interval.Query(tree, pair.R2, "any")
		}

		if len(ansR1) == 2 {
			if sameConstruct(ansR1[0].(bed.Bed).Name, ansR1[1].(bed.Bed).Name, s.DualBx) {
				twoCollapseToOne++
				ansR1 = ansR1[:1]
			}
		}
		if len(ansR2) == 2 {
			if sameConstruct(ansR2[0].(bed.Bed).Name, ansR2[1].(bed.Bed).Name, s.DualBx) {
				twoCollapseToOne++
				ansR2 = ansR2[1:]
			}
		}

		switch {
		case len(ansR1) == 1 && len(ansR2) == 1:
			if sameConstruct(ansR1[0].(bed.Bed).Name, ansR2[0].(bed.Bed).Name, s.DualBx) {
				if s.DualBx {
					nameR1 = stripDualBx(ansR1[0].(bed.Bed).Name)
				} else {
					nameR1 = ansR1[0].(bed.Bed).Name
				}
				curr = countsMap[nameR1]
				curr.rawCounts++
				if intron {
					curr.intronCount++

				}
				countsMap[nameR1] = curr
				same++
			} else {
				if s.DualBx {
					nameR1 = stripDualBx(ansR1[0].(bed.Bed).Name)
					curr = countsMap[nameR1]
					curr.rawCounts += 0.5
					curr.halfCounts++
					if intron {
						curr.intronCount += 0.5
					}
					countsMap[nameR1] = curr

					nameR2 = stripDualBx(ansR2[0].(bed.Bed).Name)
					curr = countsMap[nameR2]
					curr.rawCounts += 0.5
					curr.halfCounts++
					if intron {
						curr.intronCount += 0.5
					}
					countsMap[nameR2] = curr

					twoBxHalfCounts++
				}
			}
		case len(ansR1)+len(ansR2) == 1:
			if len(ansR1) == 1 {
				if s.DualBx {
					nameR1 = stripDualBx(ansR1[0].(bed.Bed).Name)
				} else {
					nameR1 = ansR1[0].(bed.Bed).Name
				}
				if pair.R2.QName == "" {
					oneOnlyOtherUnmap++
				} else {
					oneOnlyOtherFilterd++
				}

			} else {
				if s.DualBx {
					nameR1 = stripDualBx(ansR2[0].(bed.Bed).Name)
				} else {
					nameR1 = ansR2[0].(bed.Bed).Name
				}
				if pair.R1.RName == "" {
					oneOnlyOtherFilterd++
				}
			}

			curr = countsMap[nameR1]
			curr.rawCounts += 0.5
			curr.halfCounts++
			if intron {
				curr.intronCount += 0.5
			}
			countsMap[nameR1] = curr
			oneBxHalfCount++

		case (len(ansR1) == 1 && len(ansR2) == 2) || (len(ansR1) == 2 && len(ansR2) == 1):
			nameR1 = rescueTotalThree(ansR1, ansR2, s)
			if nameR1 != "" {
				curr = countsMap[nameR1]
				curr.rawCounts++
				if intron {
					curr.intronCount++

				}
				countsMap[nameR1] = curr
				totLen3_save++
			} else {
				totLen3_discard++
			}
		case len(ansR1) > 1 && len(ansR2) > 1:
			//fmt.Println(ansR1, ansR2)
			bothmultianddiff++
		case len(ansR1)+len(ansR2) == 0:
			none++
		default:
			//sam.WriteSamPeToFileHandle(outPe, pair)
			other++
		}
	}
	fmt.Println("Processing ", s.InSam)
	//totPair, qualIntronFilter, twoCollapseToOne, diff, same, oneOfTwo, bothmultianddiff, none, other
	fmt.Println("total read pairs: ", totPair)
	fmt.Println("reads (not pairs) filtered out by quality or intron filters: ", qualIntronFilter) //filtered out, feel good about it?
	fmt.Println("read pairs completely filtered out due to quality or intron filters: ", pairsFiltered)
	fmt.Println("reads (not pairs) that overlapped both indexes from the same construct: ", twoCollapseToOne) //counted, feel good about it
	fmt.Println("pairs in which each read maps to the same construct (opposite bx): ", same)
	fmt.Println("barcode overlaps that had >1 mismatch to ref and were filtered out: ", bxMismatch)
	fmt.Println("Read pairs in which both reads mapped to a single barcode, but to different constructs. (half counts awarded): ", twoBxHalfCounts)        // not using currently, want to see if I can rescue
	fmt.Println("Read pairs in which one read mapped to a barcode, and it's read pair didn't overlap any barcodes (half count awarded): ", oneBxHalfCount) // using currently, need to figure out why the other read doesn't map
	fmt.Println("How many of the above pairs only mapped to one total construct because the other read pair was filtered out: ", oneOnlyOtherFilterd)
	fmt.Println("How many of the above pairs only mapped to one total construct because the other read pair was unmapped by STAR: ", oneOnlyOtherUnmap)
	fmt.Println("How many of the above pairs only mapped to one total construct because the other read pair didn't overlap a Bx (intron usually?): ", oneBxHalfCount-(oneOnlyOtherFilterd+oneOnlyOtherUnmap))
	fmt.Println("Total length of three saved and discarded: ", totLen3_save, totLen3_discard)
	fmt.Println("Read pairs where both reads map to 2+ barcodes: ", bothmultianddiff) // not using need to look into
	fmt.Println("Read pairs where neither pair overlapped a barcode: ", none)         // obvi not using, haven't looked at
	fmt.Println("other, didn't match any of the other categories: ", other)

	if s.InputNorm != "" {
		inputNormalize(s.InputNorm, countsMap)
	}
	if s.ZScore != "" {
		zscoreNorm(s, countsMap)
	}
	if s.IntronStats {
		intronMath(countsMap)
	}

	return formatResults(countsMap, s)
}

func readBxHoppingFile(countsMap map[string]construct, bxHoppingFile string) {
	var cols []string
	var curr construct
	slc := fileio.Read(bxHoppingFile)
	for i := range slc {
		cols = strings.Split(slc[i], "\t")
		curr = countsMap[cols[0]]
		curr.family = cols[1]
		curr.barcode = dna.StringToBases(cols[2])
		countsMap[cols[0]] = curr
	}
}

func processSingleEndSam(s BulkOutputSeqSettings, tree map[string]*interval.IntervalNode, countsMap map[string]construct, ref []fasta.Fasta) []string {
	var ans []interval.Interval
	var zero, multi, bxMismatch, toAdd int
	var curr construct
	var intron bool
	var name string

	inSam, _ := sam.GoReadToChan(s.InSam)

	for read := range inSam {
		if read.MapQ < 1 {
			continue
		} else if !filterIntrons(read) {
			continue
		}

		intron = false
		if detectIntron(read) || detectIntron(read) {
			intron = true
		}

		if s.DualBx || s.SingleBx {
			ans, toAdd = bxOverlap(read, tree, ref, s, countsMap)
			bxMismatch += toAdd
		} else {
			ans = interval.Query(tree, read, "any")
		}

		if len(ans) == 2 {
			if sameConstruct(ans[0].(bed.Bed).Name, ans[1].(bed.Bed).Name, s.DualBx) {
				ans = ans[:1]
			}
		}

		switch len(ans) {
		case 0:
			zero++
			continue
		case 1:
			if s.DualBx {
				name = stripDualBx(ans[0].(bed.Bed).Name)
			} else {
				name = ans[0].(bed.Bed).Name
			}
			curr = countsMap[name]
			curr.rawCounts++
			if intron {
				curr.intronCount++
			}
			countsMap[name] = curr
		default:
			multi++
			continue
		}
	}

	if s.IntronStats {
		intronMath(countsMap)
	}

	if s.InputNorm != "" {
		inputNormalize(s.InputNorm, countsMap)
	}
	if s.ZScore != "" {
		zscoreNorm(s, countsMap)
	}
	return formatResults(countsMap, s)
}

func detectIntron(s sam.Sam) bool {
	for i := range s.Cigar {
		if s.Cigar[i].Op == cigar.Ns && s.Cigar[i].RunLength >= 100 {
			return true
		}
	}
	return false
}

func zscoreNorm(s BulkOutputSeqSettings, countsMap map[string]construct) {
	var ncCounts []float64
	var c construct
	nc := fileio.Read(s.ZScore)
	for i := range nc {
		ncCounts = append(ncCounts, countsMap[nc[i]].normCounts)
	}
	ncAvg := numbers.AverageFloat64(ncCounts)
	ncSD := numbers.StandardDeviationFloat64(ncCounts)
	for i := range countsMap {
		c = countsMap[i]
		c.zscore = numbers.ZScore(countsMap[i].normCounts, ncAvg, ncSD)
		countsMap[i] = c
	}

}

/*
	func bxOverlap(s sam.Sam, bedTree map[string]*interval.IntervalNode, ref []fasta.Fasta, settings BulkOutputSeqSettings, countsMap map[string]construct) ([]interval.Interval, int) {
		if s.QName == "" {
			return []interval.Interval{}, 0
		}
		var tmp, ans []interval.Interval
		var toAdd bed.Bed
		var j, bxMismatch int
		var saved string
		samBeds := convert.SamToBedExcludeIntron(s)

		for i := range samBeds {
			tmp = interval.Query(bedTree, samBeds[i], "d")
			if settings.CheckBx != "" {
				ans = append(ans, tmp...)
				continue
			}
			for j = range tmp {
				saved = checkBx(s, tmp[j].(bed.Bed), ref, countsMap)
				if saved != "" {
					toAdd = tmp[j].(bed.Bed)
					if settings.DualBx {
						toAdd.Name = saved + "_a"
					} else {
						toAdd.Name = saved
					}
					toAdd.Name = saved
					ans = append(ans, toAdd)
				} else {
					bxMismatch++
				}
			}
		}
		return ans, bxMismatch
	}
*/
func bxOverlap(s sam.Sam, bedTree map[string]*interval.IntervalNode, ref []fasta.Fasta, settings BulkOutputSeqSettings, countsMap map[string]construct) ([]interval.Interval, int) {
	if s.QName == "" {
		return []interval.Interval{}, 0
	}
	var tmp, ans []interval.Interval
	var toAdd bed.Bed
	var j, bxMismatch int
	var saved string

	tmp = interval.Query(bedTree, s, "d")
	if settings.CheckBx == "" {
		ans = append(ans, tmp...)
		return ans, 0
	}
	for j = range tmp {
		saved = checkBx(s, tmp[j].(bed.Bed), ref, countsMap)
		if saved != "" {
			toAdd = tmp[j].(bed.Bed)
			if settings.DualBx {
				toAdd.Name = saved + "_a"
			} else {
				toAdd.Name = saved
			}
			toAdd.Name = saved
			ans = append(ans, toAdd)
		} else {
			bxMismatch++
		}
	}
	return ans, bxMismatch
}

func checkBx(s sam.Sam, bxCoord bed.Bed, refs []fasta.Fasta, countsMap map[string]construct) string {
	var mismatch int
	var err error
	fmt.Println(s)
	fmt.Println(bxCoord)
	obs, err := sam.SamBedToBases(s, bxCoord)
	if err != nil {
		return ""
	}
	refbx := convert.BedToFasta([]bed.Bed{bxCoord}, refs)[0]

	for i := range obs {
		if obs[i] != dna.ToUpper(refbx.Seq[i]) {
			mismatch++
		}
		if mismatch > 1 && countsMap[s.RName].family != "" {
			return checkForChimeras(obs, countsMap, s.RName)
		}
	}
	return s.RName
}

func checkForChimeras(obs []dna.Base, countsMap map[string]construct, alnChrom string) string {
	for i := range countsMap {
		if dna.CompareSeqsIgnoreCase(countsMap[i].barcode, obs) == 0 && countsMap[alnChrom].family == countsMap[i].family {
			return countsMap[i].name
		}
	}
	return ""
}

func rescueTotalThree(ansR1, ansR2 []interval.Interval, s BulkOutputSeqSettings) string {
	var j int
	for i := range ansR1 {
		for j = range ansR2 {
			if sameConstruct(ansR1[i].(bed.Bed).Name, ansR2[j].(bed.Bed).Name, true) {
				if s.DualBx {
					return stripDualBx(ansR1[i].(bed.Bed).Name)
				} else {
					return ansR1[i].(bed.Bed).Name
				}
			}
		}
	}
	return ""
}

func stripDualBx(s string) string {
	slc := strings.Split(s, "_")
	return strings.Join(slc[:len(slc)-1], "_")
}

// sameConstruct returns true if the barcodes come from the same construct
func sameConstruct(c1, c2 string, dualBx bool) bool {
	if !dualBx {
		return c1 == c2
	}
	name1 := stripDualBx(c1)
	name2 := stripDualBx(c2)
	return name1 == name2
}

func formatResults(countsMap map[string]construct, s BulkOutputSeqSettings) []string {
	var res []string
	var base string
	var curr construct

	for i := range countsMap {
		curr = countsMap[i]
		base = fmt.Sprintf("%s\t%.1f", curr.name, curr.rawCounts)
		if s.InputNorm != "" {
			base += fmt.Sprintf("\t%.2f\t%.2f", curr.normCounts, curr.normFactor)
		}
		if s.IntronStats {
			base += fmt.Sprintf("\t%.1f\t%.2f", curr.intronCount, curr.intronPercent)
		}
		if s.ZScore != "" {
			base += fmt.Sprintf("\t%.2f", curr.zscore)
		}
		base += fmt.Sprintf("\t%d", s.RepNum)
		base += fmt.Sprintf("\t%d", curr.halfCounts)
		res = append(res, base)
	}

	return res
}

func inputNormalize(normFactor string, countsMap map[string]construct) {
	var found bool
	var counts construct
	normMap := readNormFactorFile(normFactor)
	for i := range normMap {
		counts, found = countsMap[i]
		if !found {
			counts = construct{
				name:        i,
				rawCounts:   0,
				normCounts:  0,
				normFactor:  normMap[i],
				intronCount: 0,
				zscore:      0,
			}
			countsMap[i] = construct{}
			fmt.Printf("WARNING: %s is in the input normalization file, but no counts were found. This construct will have 0 counts in the output file\n", i)
			continue
		}
		counts.normCounts = counts.rawCounts * normMap[i]
		counts.normFactor = normMap[i]
		countsMap[i] = counts
	}
	checkIfEverythingWasNormed(countsMap, normMap)
}

func checkIfEverythingWasNormed(countsMap map[string]construct, normMap map[string]float64) {
	var found bool
	var notFound []string
	for i := range countsMap {
		_, found = normMap[i]
		if !found {
			notFound = append(notFound, i)
		}
	}
	if len(notFound) != 0 {
		log.Fatalf("ERROR: These constructs have counts but were not found in the input normalization file: %v\n", notFound)
	}
}

func makeBedTreeAndMap(inBed string) (tree map[string]*interval.IntervalNode, countsMap map[string]construct) {
	var ntrvls []interval.Interval
	var c construct

	countsMap = make(map[string]construct)

	bd := bed.Read(inBed)

	for i := range bd {
		ntrvls = append(ntrvls, bd[i])
		c = construct{
			name:          bd[i].Name,
			rawCounts:     0,
			normCounts:    0,
			normFactor:    0,
			intronCount:   0,
			intronPercent: 0,
			zscore:        0,
		}
		countsMap[bd[i].Name] = c
	}

	tree = interval.BuildTree(ntrvls)

	return tree, countsMap
}

// filterIntrons returns true if they pass the intron filter
func filterIntrons(s sam.Sam) bool {
	for i := range s.Cigar {
		if s.Cigar[i].Op == 'N' && s.Cigar[i].RunLength > 700 {
			return false
		}
	}
	return true
}

func readNormFactorFile(inNorm string) map[string]float64 {
	var cols []string
	normFactorMap := make(map[string]float64)
	in := fileio.Read(inNorm)
	for i := 1; i < len(in); i++ {
		cols = strings.Split(in[i], "\t")
		normFactorMap[cols[0]] = parse.StringToFloat64(cols[3])
	}
	return normFactorMap
}

func BulkOutputSeq(s BulkOutputSeqSettings) []string {
	var ref []fasta.Fasta

	tree, countsMap := makeBedTreeAndMap(s.InBed)

	if s.DualBx {
		countsMap = collapseDualBxOutput(countsMap)
	}

	if s.CheckBx != "" {
		ref = fasta.Read(s.CheckBx)
	}

	if s.PairedEnd {
		return processPairEndSam(s, tree, countsMap, ref)
	} else {
		return processSingleEndSam(s, tree, countsMap, ref)
	}
}

func collapseDualBxOutput(countsMap map[string]construct) map[string]construct {
	var slc []string
	var c construct

	collapseMap := make(map[string]construct)

	for i := range countsMap {
		slc = strings.Split(i, "_")
		c = construct{
			name:          strings.Join(slc[:len(slc)-1], "_"),
			rawCounts:     0,
			normCounts:    0,
			normFactor:    0,
			intronCount:   0,
			intronPercent: 0,
			zscore:        0,
		}
		collapseMap[c.name] = c
	}
	if len(collapseMap)*2 != len(countsMap) {
		fmt.Printf("WARNING: Dual barcode collapsed map does not have twice as fewer constructs as the original map, some constructs may not have collapsed.")
	}
	return collapseMap
}

func intronMath(countsMap map[string]construct) {
	var curr construct
	for i := range countsMap {
		curr = countsMap[i]
		curr.intronPercent = (curr.intronCount / curr.rawCounts) * 100
		countsMap[i] = curr
	}
}
