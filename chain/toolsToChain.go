package chain

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"strings"
)

//TODO: Will move to the overlap interface once we have that set up, essentially all the functions
//if target bool is true, we select the target/refernce regions, if target bool is false we use the query
func OverlapChainBed(alpha *Chain, beta *bed.Bed, checkTarget bool) bool {
	if checkTarget {
		return targetOverlap(alpha, beta)
	} else {
		return queryOverlap(alpha, beta)
	}
}

func targetOverlap(alpha *Chain, beta *bed.Bed) bool {
	if (common.MaxInt64(int64(alpha.TStart), beta.ChromStart) < common.MinInt64(int64(alpha.TEnd), beta.ChromEnd)) && strings.Compare(alpha.TName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func queryOverlap(alpha *Chain, beta *bed.Bed) bool {
	if (common.MaxInt64(int64(alpha.QStart), beta.ChromStart) < common.MinInt64(int64(alpha.QEnd), beta.ChromEnd)) && strings.Compare(alpha.QName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func ChainToBed(ch *Chain, target bool) *bed.Bed {
	if target {
		return &bed.Bed{Chrom: ch.TName, ChromStart: int64(ch.TStart), ChromEnd: int64(ch.TEnd), Name: ch.QName, Score: int64(ch.Score)}
	} else {
		return &bed.Bed{Chrom: ch.QName, ChromStart: int64(ch.QStart), ChromEnd: int64(ch.QEnd), Name: ch.TName, Score: int64(ch.Score)}
	}
}
