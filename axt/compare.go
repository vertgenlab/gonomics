package axt

import (
	//"bufio"
	//"flag"
	//"fmt"
	"sort"
	//"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/fasta"
	//"github.com/vertgenlab/gonomics/vcf"
	//"io"
	//"io/ioutil"
	//"log"
	//"os"
	//"strconv"
	"strings"
)

func CompareRCoord(alpha *Axt, beta *Axt) int {
	if alpha.RStart < beta.RStart {
		return -1
	}
	if alpha.RStart > beta.RStart {
		return 1
	}
	if alpha.REnd < beta.REnd {
		return -1
	}
	if alpha.REnd > beta.REnd {
		return 1
	}
	return 0
}

func CompareRName(alpha *Axt, beta *Axt) int {
	return strings.Compare(alpha.RName, beta.RName)
}

func CompareQName(alpha *Axt, beta *Axt) int {
	return strings.Compare(alpha.QName, beta.QName)
}

func CompareScore(alpha *Axt, beta *Axt) int {
	if alpha.Score < beta.Score {
		return 1
	}
	if alpha.Score > beta.Score {
		return -1
	}
	return 0
}

func CompareRNameCoord(alpha *Axt, beta *Axt) int {
	compareStorage := CompareRName(alpha, beta)
	if compareStorage != 0 {
		return compareStorage
	} else {
		return CompareRCoord(alpha, beta)
	}
}

func SortByRNameCoord(axts []*Axt) {
	sort.Slice(axts, func(i, j int) bool { return CompareRNameCoord(axts[i], axts[j]) == -1 })
}

func SortByScore(axts []*Axt) {
	sort.Slice(axts, func(i, j int) bool { return CompareScore(axts[i], axts[j]) == -1 })
}
