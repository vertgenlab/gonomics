package chain

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"log"
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
	var tStart, tEnd int
	if !alpha.TStrand {
		tStart, tEnd = getSwapTCoord(alpha, true, false), getSwapTCoord(alpha, false, true)
	} else {
		tStart, tEnd = alpha.TStart, alpha.TEnd
	}
	if (common.MaxInt64(int64(tStart), beta.ChromStart) < common.MinInt64(int64(tEnd), beta.ChromEnd)) && strings.Compare(alpha.TName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func queryOverlap(alpha *Chain, beta *bed.Bed) bool {
	var qStart, qEnd int
	if !alpha.QStrand {
		qStart, qEnd = getSwapQCoord(alpha, true, false), getSwapQCoord(alpha, false, true)
	} else {
		qStart, qEnd = alpha.QStart, alpha.TEnd
	}
	if (common.MaxInt64(int64(qStart), beta.ChromStart) < common.MinInt64(int64(qEnd), beta.ChromEnd)) && strings.Compare(alpha.QName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func ChainToBed(ch *Chain, checkTarget bool) *bed.Bed {
	if checkTarget {
		return convertTargetBed(ch)
	} else {
		return convertQueryBed(ch)
	}
}

//helper functions to convert 4 different cases target, positive strand and reverse. query positive strand and reverse
func convertTargetBed(ch *Chain) *bed.Bed {
	if ch.TStrand {
		return &bed.Bed{Chrom: ch.TName, ChromStart: int64(ch.TStart), ChromEnd: int64(ch.TEnd), Name: ch.QName, Score: int64(ch.Score)}
	} else {
		return reverseTStrandBed(ch)
	}
}

func convertQueryBed(ch *Chain) *bed.Bed {
	if ch.TStrand {
		return &bed.Bed{Chrom: ch.TName, ChromStart: int64(ch.TStart), ChromEnd: int64(ch.TEnd), Name: ch.QName, Score: int64(ch.Score)}
	} else {
		return reverseQStrandBed(ch)
	}
}

//converts chain to be strictly only the negative strand
//let me know what you think about me separating these two log.Fatals
//i could do if !ch.TStrand || !ch.QStrand { but i was worry about not catching ch.TStrand && ch.QStrand case
//still runs the same way i believe
func reverseTStrandBed(ch *Chain) *bed.Bed {
	if !ch.TStrand {
		return &bed.Bed{Chrom: ch.TName, ChromStart: int64(ch.TSize - ch.TEnd), ChromEnd: int64(ch.TSize - ch.TStart), Name: ch.QName, Score: int64(ch.Score)}
	} else {
		log.Fatalf("Error: Found a target alignment with positive strand, please check input...\n")
		return nil
	}
}

//Just want to point out that it might be weird checking the false statement first
//However, the more important condition is checked first
func reverseQStrandBed(ch *Chain) *bed.Bed {
	if !ch.QStrand {
		return &bed.Bed{Chrom: ch.QName, ChromStart: int64(ch.QSize - ch.QEnd), ChromEnd: int64(ch.QSize - ch.QStart), Name: ch.TName, Score: int64(ch.Score)}
	} else {
		log.Fatalf("Error: Found a query alignment with positive strand, please check input...\n")
		return nil
	}
}
