package chain

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
)

func countTarget(ch *Chain) {
	var targetDiff int = 0
	for i := 0; i < len(ch.Alignment.Data); i++ {
		targetDiff += ch.Alignment.Data[i].TDiff + ch.Alignment.Data[i].Size

	}
	fmt.Printf("Target total: %d\n", targetDiff)
}

func totalLen(ch *Chain) {
	var total int = 0
	for i := 0; i < len(ch.Alignment.Data); i++ {
		total += ch.Alignment.Data[i].Size
	}
	fmt.Printf("Block Size: %d\n", total)
}

func printTargetQueryNum(ch *Chain) {
	t := ch.TEnd - ch.TStart
	q := ch.QEnd - ch.QStart
	fmt.Printf("Target=%d bases, Query=%d bases\n", t, q)
}

func printHeader(ch *Chain) string {
	return fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
}

func printTargetQueryBedLike(ch *Chain) string {
	return fmt.Sprintf("Target: %s %d %d Query: %s %d %d\n", ch.TName, ch.TStart, ch.TEnd, ch.QName, ch.QStart, ch.QEnd)
}
