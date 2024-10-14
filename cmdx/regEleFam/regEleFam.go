package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"strings"
)

func buildTree(axtFile string) map[string]*interval.IntervalNode {
	var intrvls []interval.Interval
	ch, _ := axt.GoReadToChan(axtFile)
	for i := range ch {
		intrvls = append(intrvls, i)
	}
	return interval.BuildTree(intrvls)
}

func main() {
	flag.Parse()

	if len(flag.Args()) != 3 {
		log.Fatalf("regEleFam selfChain.axt ocr.bed out.txt\n")

	}
	var sb strings.Builder
	var outfile *fileio.EasyWriter = fileio.EasyCreate(flag.Arg(2))
	var ans []interval.Interval
	axTree := buildTree(flag.Arg(0))
	fmt.Println("done building tree")
	regEleChan := bed.GoReadToChan(flag.Arg(1))
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
	fmt.Printf("done\n")
}
