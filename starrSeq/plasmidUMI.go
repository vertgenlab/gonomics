package starrSeq

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
)

func makeTree(inBed string) map[string]*interval.IntervalNode {
	var constructs []interval.Interval
	b := bed.Read(inBed)
	for i := range b {
		constructs = append(constructs, b[i])
	}
	return interval.BuildTree(constructs)
}

func FindUMI(s sam.Sam) string {
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
			if i+20 >= len(s.Seq) {
				return ""
			}
			return dna.BasesToString(s.Seq[i+10 : i+20])
		}
		if ans == 2 {
			if i-10 < 0 {
				return ""
			}
			bs = s.Seq[i-10 : i]
			dna.ReverseComplement(bs)
			return dna.BasesToString(bs)
		}
	}
	return ""
}

func whichConstruct(ans1, ans2 []interval.Interval) string {
	switch {
	case len(ans1) == 1, len(ans2) != 1:
		return ans1[0].(bed.Bed).Name
	case len(ans1) != 1, len(ans2) == 1:
		return ans2[0].(bed.Bed).Name
	case len(ans1) == 1, len(ans2) == 1:
		return ans1[0].(bed.Bed).Name
	}
	return ""
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
