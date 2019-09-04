package qDna

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/sam"
	"reflect"
	"runtime"
	"fmt"
	"sync"
)

type NoGapAln struct {
	Start int
	End   int
	Score float64
}

var wg sync.WaitGroup

//aligner name gsw, graph smith-waterman
// O=600 E=150
var HumanChimpTwoScoreMatrix = [][]float64{
	{90, -330, -236, -356},
	{-330, 100, -318, -236},
	{-236, -318, 100, -330},
	{-356, -236, -330, 90},
}

func QDnaScore(alpha *QBase, beta *QBase, scoreMatrix [][]float64) float64 {
	var sum float64 = 0
	a := reflect.ValueOf(alpha).Elem()
	b := reflect.ValueOf(beta).Elem()
	for x := 0; x < a.NumField(); x++ {
		for y := 0; y < b.NumField(); y++ {
			sum += scoreMatrix[x][y] * a.Field(x).Float() * b.Field(y).Float()
		}
	}
	return sum
}

func QDnaFasterScore(alpha *QBase, beta *QBase, scoreMatrix [][]float64) float64 {
	var sum float64 = 0
	//Running loop one by one to see if memory alocation speeds up calculations
	sum += float64(alpha.A) * float64(beta.A) * scoreMatrix[0][0]
	sum += float64(alpha.A) * float64(beta.C) * scoreMatrix[0][1]
	sum += float64(alpha.A)* float64(beta.G) * scoreMatrix[0][2]
	sum += float64(alpha.A) * float64(beta.T) * scoreMatrix[0][3]

	sum += float64(alpha.C) * float64(beta.A) * scoreMatrix[1][0]
	sum += float64(alpha.C) * float64(beta.C) * scoreMatrix[1][1]
	sum += float64(alpha.C) * float64(beta.G) * scoreMatrix[1][2]
	sum += float64(alpha.C) * float64(beta.T) * scoreMatrix[1][3]

	sum += float64(alpha.G) * float64(beta.A) * scoreMatrix[2][0]
	sum += float64(alpha.G) * float64(beta.C) * scoreMatrix[2][1]
	sum += float64(alpha.G) * float64(beta.G) * scoreMatrix[2][2]
	sum += float64(alpha.G) * float64(beta.T) * scoreMatrix[2][3]

	sum += float64(alpha.T) * float64(beta.A) * scoreMatrix[3][0]
	sum += float64(alpha.T) * float64(beta.C) * scoreMatrix[3][1]
	sum += float64(alpha.T) * float64(beta.G) * scoreMatrix[3][2]
	sum += float64(alpha.T) * float64(beta.T) * scoreMatrix[3][3]
	return sum
}

func PairwiseAverage(alpha *QFrag, beta *QFrag, start int64, end int64, name string) *QFrag {

	answer := &QFrag{Seq: nil, From: []*Location{&Location{Assembly: "", Chr: name, Start: start, End: end}}, Fwd: nil, Rev: nil}
	//Max or min, haven't decided if it's necessary.
	//Just trying to handle a potential error when sequences are uneven
	//length = common.Min(len(alpha.Seq), len(beta.Seq))
	for i := 0; i < len(alpha.Seq); i++ {
		tmpA := (alpha.Seq[i].A + beta.Seq[i].A) / 2
		tmpC := (alpha.Seq[i].C + beta.Seq[i].C) / 2
		tmpG := (alpha.Seq[i].G + beta.Seq[i].G) / 2
		tmpT := (alpha.Seq[i].T + beta.Seq[i].T) / 2
		answer.Seq = append(answer.Seq, &QBase{A: tmpA, C: tmpC, G: tmpG, T: tmpT})
	}
	return answer
}

func reverseCigar(alpha []align.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func UngappedAlignLen(cig []align.Cigar) int64 {
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColD {
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}

func UngappedQueryLen(cig []align.Cigar) int64 {
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColI {
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}

func GSW2(ref []*QFrag, reads []*fastq.Fastq) []*sam.SamAln {
	var answer []*sam.SamAln = make([]*sam.SamAln, len(reads))
	c := make(chan *sam.SamAln)
	runtime.GOMAXPROCS(runtime.NumCPU())

	//var wg sync.WaitGroup
	wg.Add(len(reads))
	var i, j int
	m := IndexRefSlidingWindow(ref, 20)
	for i = 0; i < len(reads); i++ {
		//go gsw(ref, reads[i], c)
		go warrior(ref, reads[i], 20, m, c)
	}
	for j = 0; j < len(reads); j++ {
		answer[j] = <-c
	}
	wg.Wait()
	return answer
}

func warrior(ref []*QFrag, read *fastq.Fastq, seed int, m map[int64][]ChrDict, c chan *sam.SamAln) {
	defer wg.Done()
	var maxScore float64 = 0
	//m := indexRef(ref, seed)
	var seedBeds []*bed.Bed
	//var revSeedBeds []*bed.Bed
	//var foundSeed bool = false
	chrSize := ChromMap(ref)
	var extension int64 = int64(maxScore/600) + int64(len(read.Seq))

	//calculate max scopre for smith-waterman
	currRead := FromFastq(read)
	//Reverse strand
	reverseFastq := fastq.ReverseComplementFastq(read)
	reverseQBase := FromFastq(reverseFastq)

	for base := 0; base < len(currRead); base++ {
		maxScore += QDnaScore(currRead[base], currRead[base], HumanChimpTwoScoreMatrix)
	}

	var bestScore float64 = 0
	var currBed bed.Bed
	var bStart, bEnd int64
	var i, j int
	for i = 0; i < len(read.Seq)-seed; i++ {
		for j = 0; j < len(m[putTogether(read.Seq[i:i+seed])]); j++ {
			//build bed slice for positions to start smith-waterman
			bStart = m[putTogether(read.Seq[i:i+seed])][j].Coord - extension
			if bStart < 0 {
				bStart = 0
			}
			bEnd = m[putTogether(read.Seq[i:i+seed])][j].Coord + extension
			if bEnd > int64(len(chrSize[m[putTogether(read.Seq[i:i+seed])][j].Chr])) {
				bEnd = int64(len(chrSize[m[putTogether(read.Seq[i:i+seed])][j].Chr]))
			}
			currBed = bed.Bed{Chrom: m[putTogether(read.Seq[i:i+seed])][j].Chr, ChromStart: bStart, ChromEnd: bEnd, Strand: true}
			seedBeds = append(seedBeds, &currBed)
			//seedBeds = append(seedBeds, &bed.Bed{Chrom: m[putTogether(read.Seq[i:i+seed])][j].Chr, ChromStart: bStart, ChromEnd: bEnd, Strand: true})
		}
	}

	bed.Sort(seedBeds)
	seedBeds = bed.MergeBeds(seedBeds)

	var score, unGappedScore float64
	var alignment []align.Cigar
	var lowRef, lowQuery, highQuery int64
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}

	for k := 0; i < len(seedBeds); k++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(chrSize[seedBeds[k].Chrom], currRead, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = seedBeds[k].Chrom
			currBest.Pos = lowRef + seedBeds[k].ChromStart

			currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
			currBest.Seq = dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery]))
			currBest.Qual = string(read.Qual[lowQuery:highQuery])
		}
	}


	
	seedBeds = nil
	var x, y int
	for x = 0; i < len(reverseFastq.Seq)-seed; x++ {
		for y = 0; y < len(m[putTogether(reverseFastq.Seq[x:x+seed])]); y++ {
			//build bed slice for positions to start smith-waterman
			bStart = m[putTogether(reverseFastq.Seq[x:x+seed])][y].Coord - extension
			if bStart < 0 {
				bStart = 0
			}
			bEnd = m[putTogether(reverseFastq.Seq[x:x+seed])][j].Coord + extension
			if bEnd > int64(len(chrSize[m[putTogether(reverseFastq.Seq[x:x+seed])][y].Chr])) {
				bEnd = int64(len(chrSize[m[putTogether(reverseFastq.Seq[x:x+seed])][y].Chr]))
			}
			currBed = bed.Bed{Chrom: m[putTogether(reverseFastq.Seq[x:x+seed])][y].Chr, ChromStart: bStart, ChromEnd: bEnd, Strand: true}
			seedBeds = append(seedBeds, &currBed)
			//seedBeds = append(seedBeds, &bed.Bed{Chrom: m[putTogether(read.Seq[i:i+seed])][j].Chr, ChromStart: bStart, ChromEnd: bEnd, Strand: true})
		}
	}
	
	seedBeds = bed.MergeBeds(seedBeds)
	bed.Sort(seedBeds)
	
	
	for z := 0; z < len(seedBeds); z++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(chrSize[seedBeds[z].Chrom], reverseQBase, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 16
			currBest.RName = seedBeds[i].Chrom
			currBest.Pos = lowRef + seedBeds[z].ChromStart

			currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
			currBest.Seq = dna.BasesToString(mostLikelySeq(reverseQBase[lowQuery:highQuery]))
			currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])
		}
	}

	//case if could not find seed for either reverse or positive stand do smithwaterman on both forward and reverse
	var unGappedMinJ, unGappedMaxJ int64
	if bestScore == 0 {
	var chr int
		for chr = 0; chr < len(ref); chr++ {
			unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[chr].Seq, currRead, HumanChimpTwoScoreMatrix)
			//fmt.Printf("%d %d\n", unGappedMinJ, unGappedMaxJ)
			extension = int64(unGappedScore / 600)
			unGappedMinJ = unGappedMinJ - extension
			if unGappedMinJ < 0 {
				unGappedMinJ = 0
			}
			unGappedMaxJ = unGappedMaxJ + extension
			if unGappedMaxJ > int64(len(ref[chr].Seq)) {
				unGappedMaxJ = int64(len(ref[chr].Seq))
			}

			//score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq[unGappedMinJ:unGappedMaxJ], currRead, HumanChimpTwoScoreMatrix, -600)
			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[chr].Seq[unGappedMinJ:unGappedMaxJ], currRead, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score
				currBest.Flag = 0
				currBest.RName = ref[chr].From[0].Chr
				currBest.Pos = lowRef + unGappedMinJ

				currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
				currBest.Seq = dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery]))
				currBest.Qual = string(read.Qual[lowQuery:highQuery])
			}
		}
		for chr = 0; chr < len(ref); chr++ {
			unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[chr].Seq, reverseQBase, HumanChimpTwoScoreMatrix)
			extension = int64(unGappedScore / 600)
			unGappedMinJ = unGappedMinJ - extension
			if unGappedMinJ < 0 {
				unGappedMinJ = 0
			}
			unGappedMaxJ = unGappedMaxJ + extension
			if unGappedMaxJ > int64(len(ref[chr].Seq)) {
				unGappedMaxJ = int64(len(ref[chr].Seq))
			}

			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[chr].Seq[unGappedMinJ:unGappedMaxJ], reverseQBase, HumanChimpTwoScoreMatrix, -600)
			//fmt.Println("Gsw v1 negative score: ", score)
			//score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[chr].Seq, reverseQBase, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score
				currBest.Flag = 16
				currBest.RName = ref[chr].From[0].Chr
				currBest.Pos = lowRef + unGappedMinJ

				currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
				currBest.Seq = dna.BasesToString(mostLikelySeq(reverseQBase[lowQuery:highQuery]))
				currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])
			}
		}
	}
	
		//seedBeds = bed.MergeBeds(seedBeds)

		//bed.Sort(seedBeds)
		//score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(chrSize[seedBeds[i].Chrom], currRead, HumanChimpTwoScoreMatrix, -600)
	
	if bestScore < 1200 {
		currBest.Flag = 4
	}

	fmt.Println(read.Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
	c <- &currBest
}

func gsw(ref []*QFrag, read *fastq.Fastq, c chan *sam.SamAln) {
	//c := make(chan *sam.SamAln)

	currRead := FromFastq(read)
	var reverseRead []*QBase
	//var query []*qDna.QBase
	//var reverse []*qDna.QBase
	var score, unGappedScore float64
	var alignment []align.Cigar
	var lowRef, lowQuery, highQuery int64
	var reverseFastq *fastq.Fastq

	var bestScore float64 = 0
	var unGappedMinJ, unGappedMaxJ, extend int64
	//var currRead []*QBase
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}
	reverseFastq = fastq.ReverseComplementFastq(read)
	reverseRead = FromFastq(reverseFastq)
	var i int
	for i = 0; i < len(ref); i++ {
		unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[i].Seq, currRead, HumanChimpTwoScoreMatrix)
		//fmt.Printf("%d %d\n", unGappedMinJ, unGappedMaxJ)
		extend = int64(unGappedScore / 600)
		unGappedMinJ = unGappedMinJ - extend
		if unGappedMinJ < 0 {
			unGappedMinJ = 0
		}
		unGappedMaxJ = unGappedMaxJ + extend
		if unGappedMaxJ > int64(len(ref[i].Seq)) {
			unGappedMaxJ = int64(len(ref[i].Seq))
		}

		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[i].Seq[unGappedMinJ:unGappedMaxJ], currRead, HumanChimpTwoScoreMatrix, -600)
		//fmt.Println("Print positive Score: ", score)
		//fmt.Println(dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery])))
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = ref[i].From[0].Chr
			currBest.Pos = lowRef + unGappedMinJ

			currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
			currBest.Seq = dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery]))
			currBest.Qual = string(read.Qual[lowQuery:highQuery])
		}

		//checking the negative strand
		unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[i].Seq, reverseRead, HumanChimpTwoScoreMatrix)
		extend = int64(unGappedScore / 600)
		unGappedMinJ = unGappedMinJ - extend
		if unGappedMinJ < 0 {
			unGappedMinJ = 0
		}
		unGappedMaxJ = unGappedMaxJ + extend
		if unGappedMaxJ > int64(len(ref[i].Seq)) {
			unGappedMaxJ = int64(len(ref[i].Seq))
		}

		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[i].Seq[unGappedMinJ:unGappedMaxJ], reverseRead, HumanChimpTwoScoreMatrix, -600)
		//fmt.Println("Print negative score: ", score)
		//fmt.Println(dna.BasesToString(mostLikelySeq(reverseRead[lowQuery:highQuery])))
		if score > bestScore {
			bestScore = score
			currBest.Flag = 16
			currBest.RName = ref[i].From[0].Chr
			currBest.Pos = lowRef + unGappedMinJ
			currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
			currBest.Seq = dna.BasesToString(reverseFastq.Seq[lowQuery:highQuery])
			currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	fmt.Println(read.Name, "\t", bestScore, "\t", currBest.Cigar)
	c <- &currBest

}

//version 1
func GSW(ref []*QFrag, reads []*fastq.Fastq) []*sam.SamAln {
	var answer []*sam.SamAln = make([]*sam.SamAln, len(reads))
	var reverseRead []*QBase
	//var query []*qDna.QBase
	//var reverse []*qDna.QBase
	var score, unGappedScore float64
	var alignment []align.Cigar
	var lowRef, lowQuery, highQuery int64
	var reverseFastq *fastq.Fastq

	var bestScore float64 = 0
	//var minI, minJ, maxJ int64
	//var qualBase []rune
	//var sequence string
	//var flag int64
	var i, j int

	var unGappedMinJ, unGappedMaxJ, extend int64
	var currRead []*QBase
	for i = 0; i < len(reads); i++ {
		currRead = FromFastq(reads[i])

		reverseFastq = fastq.ReverseComplementFastq(reads[i])
		reverseRead = FromFastq(reverseFastq)

		var currBest sam.SamAln = sam.SamAln{QName: reads[i].Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}
		bestScore = 0
		for j = 0; j < len(ref); j++ {
			//query := FromDna(reads[j].Seq, ErrorRate(reads[j].Qual))
			unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[j].Seq, currRead, HumanChimpTwoScoreMatrix)
			//fmt.Printf("%d %d\n", unGappedMinJ, unGappedMaxJ)
			extend = int64(unGappedScore / 600)
			unGappedMinJ = unGappedMinJ - extend
			if unGappedMinJ < 0 {
				unGappedMinJ = 0
			}
			unGappedMaxJ = unGappedMaxJ + extend
			if unGappedMaxJ > int64(len(ref[j].Seq)) {
				unGappedMaxJ = int64(len(ref[j].Seq))
			}

			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq[unGappedMinJ:unGappedMaxJ], currRead, HumanChimpTwoScoreMatrix, -600)
			//fmt.Println("Gsw v1 Positive score: ", score)
			//fmt.Println(dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery])))
			if score > bestScore {
				bestScore = score
				currBest.Flag = 0
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef + unGappedMinJ
				currBest.MapQ = 255
				currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
				currBest.Seq = dna.BasesToString(mostLikelySeq(currRead[lowQuery:highQuery]))
				currBest.Qual = string(reads[i].Qual[lowQuery:highQuery])
			}
			unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ = UngappedAlign(ref[j].Seq, reverseRead, HumanChimpTwoScoreMatrix)
			extend = int64(unGappedScore / 600)
			unGappedMinJ = unGappedMinJ - extend
			if unGappedMinJ < 0 {
				unGappedMinJ = 0
			}
			unGappedMaxJ = unGappedMaxJ + extend
			if unGappedMaxJ > int64(len(ref[j].Seq)) {
				unGappedMaxJ = int64(len(ref[j].Seq))
			}

			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq[unGappedMinJ:unGappedMaxJ], reverseRead, HumanChimpTwoScoreMatrix, -600)
			//fmt.Println("Gsw v1 negative score: ", score)
			//fmt.Println(dna.BasesToString(mostLikelySeq(reverseRead[lowQuery:highQuery])))
			if score > bestScore {
				bestScore = score
				currBest.Flag = 16
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef + unGappedMinJ
				currBest.MapQ = 255
				currBest.Cigar = cigar.FromString(align.PrintCigar(alignment))
				currBest.Seq = dna.BasesToString(mostLikelySeq(reverseRead[lowQuery:highQuery]))
				currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])
			}

			//answer = append(answer, &sam.SamAln{QName: qName, Flag: flag, RName: readName, Pos: minI, MapQ: mappingQ, Cigar: reds, RNext: rNext, PNext: pNext, TLen: tlen, Seq: sequence, Qual: string(qualBase), Extra: ""})
		}
		answer[i] = &currBest
	}

	return answer
}


/*
func unGapSw(alpha *QFrag, beta *fastq.Fastq, scoreMatrix [][]float64) *sam.SamAln {
	//var answer *sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, Cigar: "", RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}
	unGappedScore, unGappedMinJ, unGappedMaxJ, _, _ := UngappedAlign(alpha.Seq, beta.Seq, scoreMatrix)

	extend := int64(unGappedScore / 600)
	unGappedMinJ = unGappedMinJ - extend
	if unGappedMinJ < 0 {
		unGappedMinJ = 0
	}
	unGappedMaxJ = unGappedMaxJ + extend
	if unGappedMaxJ > int64(len(alpha.Seq)) {
		unGappedMaxJ = int64(len(alpha.Seq))
	}
	score, alignment, lowRef, highRef, lowQuery, highQuery := SmithWaterman(alpha.Seq[unGappedMinJ:unGappedMaxJ], beta.Seq, scoreMatrix, -600)
	return &sam.SamAln{QName: read.Name, Flag: 0, RName: alpha.From[0].Chr, Pos: lowRef + unGappedMinJ, MapQ: 255, Cigar: align.PrintCigar(alignment), RNext: "*", PNext: 0, TLen: 0, Seq: dna.BasesToString(beta.Seq), Qual: beta.Qual, Extra: fmt.Sprintf("%.2f", score)}
}*/

/*
func GSW(ref []*QFrag, reads []*fastq.Fastq) []*sam.SamAln {
	var answer []*sam.SamAln = make([]*sam.SamAln, len(reads))
	var reverseRead []*QBase
	//var query []*qDna.QBase
	//var reverse []*qDna.QBase
	var score float64
	var alignment []align.Cigar
	var lowRef, lowQuery, highQuery int64
	var reverseFastq *fastq.Fastq

	var bestScore float64 = 0
	//var minI, minJ, maxJ int64
	//var qualBase []rune
	//var sequence string
	//var flag int64
	var i, j int
	var currRead []*QBase
	//var route =  make([]align.Cigar, 1)
	for i = 0; i < len(reads); i++ {
		currRead = FromFastq(reads[i])

		reverseFastq = fastq.ReverseComplementFastq(reads[i])
		reverseRead = FromFastq(reverseFastq)

		var currBest sam.SamAln = sam.SamAln{QName: reads[i].Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, Cigar: "", RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}
		bestScore = 0
		for j = 0; j < len(ref); j++ {

			//query := FromDna(reads[j].Seq, ErrorRate(reads[j].Qual))

			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq, currRead, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score

				currBest.Flag = 0
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef
				currBest.MapQ = 255
				currBest.Cigar = align.PrintCigar(alignment)
				currBest.Seq = dna.BasesToString(reads[i].Seq[lowQuery:highQuery])
				currBest.Qual = string(reads[i].Qual[lowQuery:highQuery])

			}
			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq, reverseRead, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score

				currBest.Flag = 16
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef
				currBest.MapQ = 255
				currBest.Cigar = align.PrintCigar(alignment)
				currBest.Seq = dna.BasesToString(reverseFastq.Seq[lowQuery:highQuery])
				currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])

			}
			answer[i] = &currBest
			//answer = append(answer, &sam.SamAln{QName: qName, Flag: flag, RName: readName, Pos: minI, MapQ: mappingQ, Cigar: reds, RNext: rNext, PNext: pNext, TLen: tlen, Seq: sequence, Qual: string(qualBase), Extra: ""})
		}
	}
	return answer
}*/

/*
func UngappedAlign(alpha []*QBase, beta []*QBase, alphaOffset int, betaOffset int, scoreMatrix [][]float64) float64 {
}*/
