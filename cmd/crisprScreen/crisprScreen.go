package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

type guide struct {
	name   string
	seq    string
	counts int
	cells  [][]dna.Base
	proper bool
}

func checkCell(g guide, peFq fastq.SingleCellPair) guide {
	var found bool = false
	for i := range g.cells {
		if !dna.SeqsAreSimilar(peFq.Bx, g.cells[i], 2) {
			continue
		} else {
			found = true
			break
		}
	}
	if !found {
		g.cells = append(g.cells, peFq.Bx)
	}
	return g
}

func readGuides(guidesFile string) (map[string]guide, []dna.Base) {
	var cols []string
	var g guide
	mp := make(map[string]guide)
	seqs := fileio.Read(guidesFile)
	for i := 1; i < len(seqs); i++ {
		cols = strings.Split(seqs[i], ",")
		g = guide{name: cols[0], seq: cols[4], counts: 0, cells: [][]dna.Base{}, proper: true}
		mp[g.seq] = g
	}
	return mp, dna.StringToBases(cols[3][0:20])
}

func analyze(mp map[string]guide, outDir string) {
	var j int
	var bx string
	cellOut := fileio.EasyCreate(outDir + "guidesPerCell.txt")
	guideOut := fileio.EasyCreate(outDir + "guideCount.txt")
	cellMap := make(map[string][]string)
	for i := range mp {
		fileio.WriteToFileHandle(guideOut, fmt.Sprintf("%s\t%s\t%d", mp[i].name, i, mp[i].counts))
		for j = range mp[i].cells {
			bx = dna.BasesToString(mp[i].cells[j])
			cellMap[bx] = append(cellMap[bx], mp[i].name)
		}
	}
	err := guideOut.Close()
	exception.PanicOnErr(err)
	for i := range cellMap {
		fileio.WriteToFileHandle(cellOut, fmt.Sprintf("%s\t%d\t%s", i, len(cellMap[i]), strings.Join(cellMap[i], ",")))
	}
	err = cellOut.Close()
	exception.PanicOnErr(err)
}

func charGuides(r1, r2, guidesFile, outDir string) {
	var noProto, c, u, p int
	var seq string
	var found bool
	var g guide

	mp, recogSeq := readGuides(guidesFile)
	fqChan := make(chan fastq.SingleCellPair)
	go fastq.ReadToChanSingleCellPair(r1, r2, 16, 10, fqChan)

	for i := range fqChan {
		c++
		if !dna.SeqsAreSimilar(recogSeq, i.Reads.Rev.Seq[34:54], 2) {
			noProto++
			continue
		}
		seq = dna.BasesToString(i.Reads.Rev.Seq[54:74])
		g, found = mp[seq]
		if found {
			g = checkCell(g, i)
			g.counts++
			mp[seq] = g
			if g.proper {
				p++
			}
		} else {
			g = guide{name: fmt.Sprintf("unknown%d", u), seq: seq, cells: [][]dna.Base{i.Bx}, counts: 0, proper: false}
			mp[seq] = g
			u++
		}
	}
	fmt.Println(noProto)
	fmt.Println("total reads assessed: ", c)
	fmt.Printf("Protospacer recognized in %f percent of cells\n", float64(c-noProto)/float64(c)*100)
	fmt.Println("recognized guide seq: ", p)
	analyze(mp, outDir)
}

func analyzeBamMap(mp map[string][]string) {
	cellOut := fileio.EasyCreate("guidesPerCell.txt")
	for i := range mp {
		fileio.WriteToFileHandle(cellOut, fmt.Sprintf("%s\t%d\t%s", i, len(mp[i]), strings.Join(mp[i], ",")))
	}
	err := cellOut.Close()
	exception.PanicOnErr(err)
}

func getBx(s sam.Sam) string {
	fields := strings.Split(s.QName, ":")
	bxUmi := fields[7]
	fmt.Println(bxUmi[0:16])
	return bxUmi[0:16]
}

func updateMap(mp map[string][]string, bx string, gRNA string) {
	var redundant bool = false
	_, found := mp[bx]
	if !found {
		mp[bx] = []string{gRNA}
		return
	} else {
		for i := range mp[bx] {
			if mp[bx][i] == gRNA {
				redundant = true
				break
			}
		}
	}
	if redundant == false {
		mp[bx] = append(mp[bx], gRNA)
	}
}

func guidesFromSam(bamFile, bedFile string) {
	var ans []sam.Sam
	var j int
	var bx string

	beds := bed.Read(bedFile)
	bamRead, _ := sam.OpenBam(bamFile)
	bai := sam.ReadBai(bamFile + ".bai")

	mp := make(map[string][]string)

	for i := range beds {
		ans = sam.SeekBamRegion(bamRead, bai, beds[i].Chrom, uint32(beds[i].ChromStart), uint32(beds[i].ChromEnd))
		for j = range ans {
			bx = getBx(ans[j])
			updateMap(mp, bx, beds[i].Name)
		}
	}
	analyzeBamMap(mp)
}

func main() {
	flag.Parse()
	//r1 := flag.Arg(0)
	//r2 := flag.Arg(1)
	bamFile := flag.Arg(0)
	bedFile := flag.Arg(1)
	//guidesFile := flag.Arg(2)
	//outDir := flag.Arg(3)

	//charGuides(r1, r2, guidesFile, outDir)
	guidesFromSam(bamFile, bedFile)
}
