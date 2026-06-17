package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"sort"
	"strings"
)

func PlasmidUMI(inSam, inBed, outFile string) {
	var constructR1, constructR2 []interval.Interval
	var construct, umi string
	var found bool
	mp := make(map[string]map[string]int)
	tree := makeTree(inBed)
	ch, _ := sam.GoReadSamPeToChan(inSam)
	for pe := range ch {
		constructR1 = interval.Query(tree, pe.R1, "any")
		constructR2 = interval.Query(tree, pe.R2, "any")
		construct = whichConstruct(constructR1, constructR2)
		umi = findUMI(pe.R1)
		if umi == "" {
			umi = findUMI(pe.R2)
		}
		if umi == "" {
			continue
		}
		_, found = mp[construct]
		if !found {
			mp[construct] = make(map[string]int)
		}
		mp[construct][umi]++
	}
	writeMap(outFile, mp)
}

func writeMap(outFile string, mp map[string]map[string]int) {
	var slc []string
	o := fileio.EasyCreate(outFile)
	fileio.WriteToFileHandle(o, "construct\tumi\tcount")
	for i := range mp {
		for j := range mp[i] {
			slc = append(slc, fmt.Sprintf("%s\t%s\t%d", i, j, mp[i][j]))
		}
	}
	var a, b []string
	sort.Slice(slc, func(i, j int) bool {
		a = strings.Fields(slc[i])
		b = strings.Fields(slc[j])
		if a[0] > b[0] {
			return false
		} else if a[0] < b[0] {
			return true
		} else {
			return parse.StringToInt(a[2]) > parse.StringToInt(b[2])
		}
	})

	for i := range slc {
		fileio.WriteToFileHandle(o, slc[i])
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func makeTree(inBed string) map[string]*interval.IntervalNode {
	var constructs []interval.Interval
	b := bed.Read(inBed)
	for i := range b {
		constructs = append(constructs, b[i])
	}
	return interval.BuildTree(constructs)
}

func findUMI(s sam.Sam) string {
	var test, find, bs []dna.Base
	var ans int

	find = dna.StringToBases("GCATGCGGAT")

	for i := 0; i < len(s.Seq)-9; i++ {
		test = s.Seq[i : i+10]
		ans = dna.CompareSeqsWithRevComp(test, find)
		if ans == -1 {
			continue
		}
		if ans == 1 {
			if i+16 >= len(s.Seq) {
				return ""
			}
			return dna.BasesToString(s.Seq[i+10 : i+16])
		}
		if ans == 2 {
			if i-6 < 0 {
				return ""
			}
			bs = s.Seq[i-6 : i]
			dna.ReverseComplement(bs)
			return dna.BasesToString(bs)
		}
	}
	return ""
}

func whichConstruct(ans1, ans2 []interval.Interval) string {
	switch {
	case len(ans1) == 1 && len(ans2) != 1:
		return ans1[0].(bed.Bed).Name
	case len(ans1) != 1 && len(ans2) == 1:
		return ans2[0].(bed.Bed).Name
	case len(ans1) == 1 && len(ans2) == 1:
		return ans1[0].(bed.Bed).Name
	}
	return "unclear"
}

/*
func CountPlasmidUMI(inSam, inBed, outFile string) {
	var ans1, ans2 []interval.Interval
	var curr, poss1, poss2 string

	mp := make(map[string]map[string]int)

	ch, _ := sam.GoReadSamPeToChan(inSam)
	constructs := makeTree(inBed)

	for i := range ch {
		ans1 = interval.Query(constructs, i.A, "any")
		ans1 = interval.Query(constructs, i.B, "any")
		if len(ans1) == 1 || len(ans2) == 1 {
			poss1 = findUMI(i.A)
			poss2 = findUMI(i.B)
			if (poss1 != "" && poss2 != "") || (poss1 == "" && poss2 == "") {
				continue
			}
			curr = whichConstruct(ans1, ans2)
			if poss1 != "" {
				mp[curr][poss1] += 1
			}
			if poss2 != "" {
				mp[curr][poss2] += 1
			}
		}
	}
	for i := range mp {
		fmt.Println(i)
		for j := range mp[i] {
			fmt.Println(j, mp[i][j])
		}
	}
}*/
