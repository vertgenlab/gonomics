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

func main() {
	var inBed, liftBed []bed.Bed
	var axtOverlap []interval.Interval
	var lifted bed.Bed = bed.Bed{Score: 0, FieldsInitialized: 5}

	flag.Parse()
	if len(flag.Args()) != 5 {
		log.Fatalf("~/go/bin/liftWithAxt in.axt ocr.bed chrom.Sizes out.lift.bed out.liftANDmerge.bed")
	}
	chromSizes := chromInfo.ReadToMap(flag.Arg(2))
	fmt.Println("built chrom sizes tree")
	axtTree := buildAxTree(flag.Arg(0))
	fmt.Println("build axt tree")
	in := bed.GoReadToChan(flag.Arg(1))
	fmt.Println("reading bed")
	for rec := range in {
		inBed = append(inBed, rec)
		axtOverlap = interval.Query(axtTree, rec, "di")
		for i := range axtOverlap {
			lifted.Chrom, lifted.ChromStart, lifted.ChromEnd = lift.LiftCoordinatesWithAxt(axtOverlap[i].(axt.Axt), rec, chromSizes[axtOverlap[i].(axt.Axt).QName].Size)
			lifted.Name = rec.Name + fmt.Sprintf("_lift%d", i)
			liftBed = append(liftBed, lifted)
		}
	}
	fmt.Println("removing self overlaps")
	outBed := removeSelfOverlaps(inBed, liftBed)
	fmt.Println("merging and writing")
	bed.Write(flag.Arg(3), outBed)
	mergedBed := bed.MergeBedsKeepNames(outBed)
	fmt.Println("length of mergedBed: ", len(mergedBed))
	fmt.Println("Re-lifting homologous")
	mergedBed = reLiftHomologous(mergedBed, axtTree, chromSizes)
	fmt.Println("length of mergedBed: ", len(mergedBed))
	bed.Write(flag.Arg(4), mergedBed)
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
			}
		}
		if ocr {
			continue
		}
		names = append([]string{fmt.Sprintf("homologousElement_%d", c)}, names...)
		c++
		mergedBed[i].Name = strings.Join(names, ",")
		homologous = append(homologous, mergedBed[i])
		bedMap[names[0]] = mergedBed[i]
	}

	fmt.Println("adding names to bed map")
	bedTree := buildBedTree(mergedBed)
	fmt.Println("build bed tree")
	var axtForLift, existingNodes []interval.Interval
	fmt.Println(len(homologous))
	for i := range homologous {
		fmt.Println(i)
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
				bd.Name = fmt.Sprintf("%s,%s_lift%d", bd.Name, names[0], c)
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
