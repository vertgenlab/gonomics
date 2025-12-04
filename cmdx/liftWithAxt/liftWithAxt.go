package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"log"
	"strings"
)

func buildAxTree(file string) map[string]*interval.IntervalNode {
	var intrvls []interval.Interval
	ch, _ := axt.GoReadToChan(file)
	for i := range ch {
		intrvls = append(intrvls, i)
	}
	return interval.BuildTree(intrvls)
}

func changeDelim(mergedBed []bed.Bed) {
	for i := range mergedBed {
		mergedBed[i].Annotation = []string{strings.Join(mergedBed[i].Annotation, ",")}
	}
}

func reLiftHomologous(mergedBed []bed.Bed, axtTree map[string]*interval.IntervalNode, chromSizes map[string]chromInfo.ChromInfo) []bed.Bed {
	var j, c, k, d int
	var names []string
	var homologous []bed.Bed
	var lifted, bd bed.Bed
	var ocr, pass, found bool
	var name string

	bedMap := make(map[string]bed.Bed)

	for i := range mergedBed {
		ocr = false
		names = strings.Split(mergedBed[i].Name, ",")
		for j = range names {
			if !strings.Contains(names[j], "lift") {
				bedMap[names[j]] = mergedBed[i]
				ocr = true
				break
			}
		}
		if ocr {
			continue
		}
		names = append([]string{fmt.Sprintf("homologousElement_%d", c)}, names...)
		mergedBed[i].Annotation = append([]string{"0.00"}, mergedBed[i].Annotation...)
		c++
		mergedBed[i].Name = strings.Join(names, ",")
		homologous = append(homologous, mergedBed[i])
		bedMap[names[0]] = mergedBed[i]
	}

	fmt.Println("adding names to bed map")
	bedTree := buildBedTree(mergedBed)
	fmt.Println("build bed tree")
	var axtForLift, existingNodes []interval.Interval
	for i := range homologous {
		names = strings.Split(homologous[i].Name, ",")
		c = 0
		axtForLift = interval.Query(axtTree, homologous[i], "di")
		for j = range axtForLift {
			lifted.Chrom, lifted.ChromStart, lifted.ChromEnd = lift.LiftCoordinatesWithAxt(axtForLift[j].(axt.Axt), homologous[i], chromSizes[homologous[j].Chrom].Size)
			existingNodes = interval.Query(bedTree, lifted, "any")
			if len(existingNodes) == 0 {
				d++
				continue
			}
			for k = range existingNodes {
				pass, name = overlapType(existingNodes[k].(bed.Bed))
				if !pass {
					continue
				}
				bd, found = bedMap[name]
				if !found {
					log.Fatalf("we should have found this node: %s", name)
				}
				if lift.AxtPercentIdentityInInterval(axtForLift[j].(axt.Axt), homologous[i]) < 70 {
					continue
				}
				bd.Name = fmt.Sprintf("%s,%s_lift%d", bd.Name, names[0], c)
				bd.Annotation = append(bd.Annotation, fmt.Sprintf("%.2f", lift.AxtPercentIdentityInInterval(axtForLift[j].(axt.Axt), homologous[i])))
				bedMap[name] = bd
				c++
			}
		}
	}
	fmt.Println("This many new nodes would be added: ", d)
	return bedMapToSlice(bedMap)
}

func bedMapToSlice(bedMap map[string]bed.Bed) []bed.Bed {
	var bedSlice []bed.Bed
	for i := range bedMap {
		bedSlice = append(bedSlice, bedMap[i])
	}
	bed.SortByCoord(bedSlice)
	return bedSlice
}

func overlapType(overlap bed.Bed) (bool, string) {
	names := strings.Split(overlap.Name, ",")
	for i := range names {
		if strings.Contains(names[i], "homologous") {
			return true, names[i]
		}
	}
	return false, ""
}

func removeSelfOverlaps(inBed, liftBed []bed.Bed) []bed.Bed {
	var ans []interval.Interval
	var j int
	var pass bool
	inBedTree := buildBedTree(inBed)
	for i := range liftBed {
		pass = true
		ans = interval.Query(inBedTree, liftBed[i], "any")
		if len(ans) == 0 {
			inBed = append(inBed, liftBed[i])
			continue
		}
		for j = range ans {
			if strings.Contains(liftBed[i].Name, ans[j].(bed.Bed).Name+"_") {
				pass = false
				continue
			}
		}
		if pass {
			inBed = append(inBed, liftBed[i])
		}
	}
	return inBed
}

func buildBedTree(inBed []bed.Bed) map[string]*interval.IntervalNode {
	var intervals []interval.Interval
	for i := range inBed {
		intervals = append(intervals, inBed[i])
	}
	return interval.BuildTree(intervals)
}

func usage() {
	fmt.Print("liftWithAxt -- Intended for use in a noncoding element clustering analysis. " +
		"This program takes in an AXT alignment file (usually a whole-genome self alignment), a set of regions in bed format the user" +
		"wants to lifted to paralogous regions, and a chrom sizes file. Bed regions will be lifted if:\n" +
		"\t1. The bed region is completely contained within an axt record\n" +
		"\t2. The bed region has >= 70% identity to the lift region (both base mismatches and gapped bases count as penalties)\n" +
		"This program does not currently support partial lifts. A bed file with all original regions and lifts will produced and a second file " +
		"will be produced after merging overlapping regions. Lifted regions that don't overlap any original regions will be re-lifted, to draw additional edges for graph-based clustering.\n" +
		"Usage:\n" +
		"liftWithAxt in.axt elements.bed chrom.Sizes out.lift.bed out.liftAndMerge.bed\n")
	flag.PrintDefaults()
}

func main() {
	var inBed, liftBed []bed.Bed
	var axtOverlap []interval.Interval
	var lifted bed.Bed = bed.Bed{Score: 0, Strand: '.', FieldsInitialized: 7}

	flag.Parse()
	flag.Usage = usage
	if len(flag.Args()) != 5 {
		usage()
		log.Fatalf("ERROR: expected 5 arguments but got %d\n", len(flag.Args()))
	}

	chromSizes := chromInfo.ReadToMap(flag.Arg(2))
	fmt.Println("built chrom sizes tree")
	axtTree := buildAxTree(flag.Arg(0))
	fmt.Println("build axt tree")

	in := bed.GoReadToChan(flag.Arg(1))
	fmt.Println("reading bed")
	for rec := range in {
		rec.Annotation = []string{"0.00"}
		rec.FieldsInitialized = 7
		rec.Strand = '.'
		inBed = append(inBed, rec)
		axtOverlap = interval.Query(axtTree, rec, "di")
		for i := range axtOverlap {
			if lift.AxtPercentIdentityInInterval(axtOverlap[i].(axt.Axt), rec) < 70 {
				continue
			}
			lifted.Chrom, lifted.ChromStart, lifted.ChromEnd = lift.LiftCoordinatesWithAxt(axtOverlap[i].(axt.Axt), rec, chromSizes[axtOverlap[i].(axt.Axt).QName].Size)
			lifted.Name = rec.Name + fmt.Sprintf("_lift%d", i)
			lifted.Annotation = []string{fmt.Sprintf("%.2f", lift.AxtPercentIdentityInInterval(axtOverlap[i].(axt.Axt), rec))}
			liftBed = append(liftBed, lifted)
		}
	}

	fmt.Println("removing self overlaps")
	outBed := removeSelfOverlaps(inBed, liftBed)
	fmt.Println("writing lifted bed")
	bed.Write(flag.Arg(3), outBed)
	fmt.Println("merging bed")
	mergedBed := bed.MergeBedsKeepNamesAndAnnotations(outBed)
	fmt.Println("length of mergedBed: ", len(mergedBed))
	//bed.Write("families/mergedBedPreReLift.bed", mergedBed)
	fmt.Println("Re-lifting homologous")
	mergedBed = reLiftHomologous(mergedBed, axtTree, chromSizes)
	fmt.Println("length of mergedBed: ", len(mergedBed))
	changeDelim(mergedBed)
	bed.Write(flag.Arg(4), mergedBed)
}
