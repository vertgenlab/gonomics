package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"strings"
)

func buildTree(axtFile string) map[string]*interval.IntervalNode {
	var intrvls []interval.Interval
	var c int
	ch, _ := axt.GoReadToChan(axtFile)
	for i := range ch {
		c++
		fmt.Println(c)
		intrvls = append(intrvls, i)
	}
	return interval.BuildTree(intrvls)
}

func main() {
	var sb strings.Builder
	var dir string = "/Users/sethweaver/Downloads/hsSD/"
	var outfile *fileio.EasyWriter = fileio.EasyCreate(dir + "regEleFam.txt")
	var ans []interval.Interval
	axTree := buildTree(dir + "hs1.selfChain.axt")
	fmt.Println("done building tree")
	regEleChan := bed.GoReadToChan(dir + "h9_atac_sharedPeaks.bed")
	for i := range regEleChan {
		sb.Reset()
		sb.WriteString(bed.ToString(i, 5))
		ans = interval.Query(axTree, i, "any")
		sb.WriteString("\t")
		for j := range ans {
			sb.WriteString(ans[j].(axt.Axt).QName + ":")
			sb.WriteString(fileio.IntToString(ans[j].(axt.Axt).QStart) + "-")
			sb.WriteString(fileio.IntToString(ans[j].(axt.Axt).QEnd) + "_")
			sb.WriteString(fileio.IntToString(ans[j].(axt.Axt).Score) + ";")
		}
		fileio.WriteToFileHandle(outfile, sb.String())
	}
}
