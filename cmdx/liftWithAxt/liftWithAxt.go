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
	bed.Write(flag.Arg(4), mergedBed)
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
			fmt.Println(ans[j].(bed.Bed).Name+"_", liftBed[i].Name)
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
