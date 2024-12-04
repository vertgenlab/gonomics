package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"github.com/vertgenlab/gonomics/numbers/parse"
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

func buildBedTree(file string) (map[string]*interval.IntervalNode, []interval.Interval) {
	var intrvls []interval.Interval
	ch := bed.GoReadToChan(file)
	for i := range ch {
		intrvls = append(intrvls, i)
	}
	return interval.BuildTree(intrvls), intrvls
}

func updateBedTree(intervalSlice []interval.Interval, new interval.Interval) (map[string]*interval.IntervalNode, []interval.Interval) {
	fmt.Println("in updateBedTree")
	newSlice := append(intervalSlice, new)
	return interval.BuildTree(newSlice), newSlice
}

func main() {
	flag.Parse()

	if len(flag.Args()) != 5 {
		log.Fatalf("regEleFam proportion selfChain.axt ocr.bed chrom.sizes out.txt\n")
	}

	var sb strings.Builder
	var chromSizes map[string]chromInfo.ChromInfo = chromInfo.ReadToMap(flag.Arg(3))
	var outfile *fileio.EasyWriter = fileio.EasyCreate(flag.Arg(4))
	var chainOverlap, bedOverlap, pass []interval.Interval
	var lifted bed.Bed = bed.Bed{FieldsInitialized: 3}
	var c int
	axTree := buildAxTree(flag.Arg(1))
	bedTree, bedIntervals := buildBedTree(flag.Arg(2))
	fmt.Println("done building trees")
	regEleChan := bed.GoReadToChan(flag.Arg(2))
	fmt.Println("bed chan read")
	for i := range regEleChan {
		sb.Reset()
		sb.WriteString(bed.ToString(i, 5))
		chainOverlap = interval.Query(axTree, i, "di")
		sb.WriteString("\t")
		c = 0
		for j := range chainOverlap {
			lifted.Chrom, lifted.ChromStart, lifted.ChromEnd = lift.LiftCoordinatesWithAxt(chainOverlap[j].(axt.Axt), i, chromSizes[chainOverlap[j].(axt.Axt).QName].Size)
			bedOverlap = interval.Query(bedTree, lifted, "any")
			pass = []interval.Interval{}
			for k := range bedOverlap {
				if interval.OverlapProportionRecursive(lifted, bedOverlap[k], parse.StringToFloat64(flag.Arg(0))) {
					pass = append(pass, bedOverlap[k])
				}
			}
			switch len(pass) {
			case 0:
				fmt.Println("didn't overlap anything...shouldn't happen")
				lifted.Name = i.Name + fmt.Sprintf("_lift%d", c)
				c++
				bedTree, bedIntervals = updateBedTree(bedIntervals, lifted)
				sb.WriteString(lifted.Chrom + "_" + fileio.IntToString(lifted.ChromStart) + "_" + fileio.IntToString(lifted.ChromEnd) + "_" + lifted.Name + ";")
			case 1:
				for l := range pass {
					sb.WriteString(pass[l].(bed.Bed).Chrom + "_" + fileio.IntToString(pass[l].(bed.Bed).ChromStart) + "_" + fileio.IntToString(pass[l].(bed.Bed).ChromEnd) + "_" + pass[l].(bed.Bed).Name + ";")
				}
			}

		}
		fileio.WriteToFileHandle(outfile, sb.String())
	}
	exception.PanicOnErr(outfile.Close())
	fmt.Printf("done\n")
}
