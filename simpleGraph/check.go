package simpleGraph

import (
	"strings"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/sam"
	"log"

)

func checkAlignment(simSam *sam.SamAln) bool {
	var answer bool = false
	words := strings.Split(simSam.QName, "_")
	var ref int64
	//var query int64
	if common.StringToInt64(words[1]) == simSam.Pos {
		ref, _ = cigar.ParseCigarAll(simSam.Cigar)
		if simSam.Pos + ref == common.StringToInt64(words[2]) {
			answer = true
		}
	}
	return answer
}

func CheckAnswers(query []*sam.SamAln) {
	var yes, no int64 = 0, 0
	for i := 0; i < len(query);i++ {
		if checkAlignment(query[i]) {
			yes++
			//log.Printf(sam.SamAlnToString(query[i]))
		} else {
			no++
			//log.Printf("This did not map:\n%s\n", sam.SamAlnToString(query[i]))
		}
	}
	log.Printf("Total number of reads aligned: %d...\n\tNumber of reads correctly aligned: %d...\n\tNumber of reads mismapped: %d...\n", len(query), yes, no)

}