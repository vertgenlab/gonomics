package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sort"
)

type settings struct {
	inSam    string
	inBed    string
	outFile  string
	norm     bool
	multiple bool
}

func createSliceInterval(s settings) ([]interval.Interval, float64, map[string]float64) {
	var slc []interval.Interval
	b := bed.Read(s.inBed)
	total := bed.TotalSize(b)
	mp := make(map[string]float64)
	for i := range b {
		mp[b[i].Name] = 0
		slc = append(slc, b[i])
	}
	return slc, float64(total) / float64(len(b)), mp
}

func createTree(b []interval.Interval) map[string]*interval.IntervalNode {
	return interval.BuildTree(b)
}

func writeMap(mp map[string]float64, avg float64, s settings) {
	var slc []string
	if s.norm {
		normalize(mp, avg, s)
	}
	for i := range mp {
		slc = append(slc, fmt.Sprintf("%s\t%f", i, mp[i]))
	}
	sort.Strings(slc)
	slc = append([]string{"bedRegion\tcount"}, slc...)
	fileio.Write(s.outFile, slc)
}

func normalize(mp map[string]float64, avg float64, s settings) {
	var sz int
	b := bed.Read(s.inBed)
	for i := range b {
		sz = size(b[i])
		mp[b[i].Name] = mp[b[i].Name] * avg / float64(sz)
	}
}

func size(b bed.Bed) int {
	return numbers.AbsInt(b.ChromStart - b.ChromEnd)
}

func bedCountSam(s settings) {
	var ans []interval.Interval
	slc, avg, mp := createSliceInterval(s)
	tree := createTree(slc)

	ch, _ := sam.GoReadToChan(s.inSam)
	for i := range ch {
		ans = interval.Query(tree, i, "any")
		switch len(ans) {
		case 0:
			continue
		case 1:
			mp[ans[0].(bed.Bed).Name]++
		default:
			if !s.multiple {
				continue
			}
			for j := range ans {
				mp[ans[j].(bed.Bed).Name]++
			}
		}
	}
	writeMap(mp, avg, s)
}

func usage() {

}

func main() {
	var norm *bool = flag.Bool("norm", false, "Normalize counts to the average size of the bed elements")
	var multipleAssign *bool = flag.Bool("multipleAssign", false, "Allow one sam record contribute to the counts of"+
		"multiple bed intervals. Default is multiple asignment reads are discarded")
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Expected %d arguments, got %d", expectedNumArgs, len(flag.Args()))
	}
	var s settings = settings{
		inSam:    flag.Arg(0),
		inBed:    flag.Arg(1),
		outFile:  flag.Arg(2),
		norm:     *norm,
		multiple: *multipleAssign,
	}
	bedCountSam(s)
}
