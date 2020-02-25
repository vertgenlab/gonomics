package simpleGraph

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"os"
	"strings"
	//"log"
)

func GraphSmithWaterman(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, bestScore int64
	var leftPath, rightPath []uint32

	var currScore int64 = 0
	perfectScore := perfectMatch(read)
	extension := int(perfectScore/600) + len(read.Seq)

	var currRead *fastq.Fastq = nil
	var seeds []*SeedDev = findSeedsInMapDev(seedHash, read, seedLen, stepSize, true)
	seeds = GraphDictionary(seeds, gg, read)

	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInMapDev(seedHash, revCompRead, seedLen, stepSize, false)
	revCompSeeds = GraphDictionary(revCompSeeds, gg, revCompRead)
	seeds = append(seeds, revCompSeeds...)
	//CompareBlastScore(seeds, read)
	var tailSeed *SeedDev
	//, m [][]int64, trace [][]rune
	//m, trace := SwMatrixSetup(int64(extension+1)
	var seedScore int64
	for i = 0; i < len(seeds) && isSeedBetter(i, seeds, bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		//log.Printf("seed hit: %d, len=%d\n", i, sumLen(seeds[i]))
		tailSeed = toTail(seeds[i])
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}

		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, extension, currRead.Seq[:seeds[i].QueryStart], m, trace)
		//log.Printf("NodeLen=%d, TargetStart=%d, length=%d\n", len(gg.Nodes[tailSeed.TargetId].Seq), tailSeed.TargetStart, tailSeed.Length)
		seedScore = BlastSeed(seeds[i], currRead)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension, currRead.Seq[tailSeed.QueryStart+tailSeed.Length:], m, trace)

		currScore = leftScore + seedScore + rightScore

		if currScore > bestScore {
			//log.Printf("Index: %d left=%d, seed=%d, right=%d\n", i, leftScore, seedScore, rightScore)
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.RName = fmt.Sprintf("%s_%d", gg.Nodes[seeds[i].TargetId].Name, seeds[i].TargetId)
			currBest.Pos = int64(minTarget) + 1
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(seeds[i])), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(minQuery, len(currRead.Seq), currBest.Cigar)
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(seeds[i])), rightPath), gg)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	return &currBest
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

func AddSClip(front int, lengthOfRead int, cig []*cigar.Cigar) []*cigar.Cigar {
	var runLen int64 = cigar.QueryLength(cig)
	if runLen < int64(lengthOfRead) {
		answer := make([]*cigar.Cigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, &cigar.Cigar{RunLength: int64(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+int(cigar.QueryLength(cig)) < lengthOfRead {
			answer = append(answer, &cigar.Cigar{RunLength: int64(lengthOfRead-front) - runLen, Op: 'S'})
		}
		return answer
	} else {
		return cig
	}

}

//perfect match
func perfectMatch(read *fastq.Fastq) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func scoreSeed(seed *SeedDev, read *fastq.Fastq) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func BasicHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	var words string = "@HD\tVN:1.6\tSO:unsorted\n"
	for i := 0; i < len(ref); i++ {
		words += fmt.Sprintf("@SQ\tSN:%s_%d\tLN:%d\n", ref[i].Name, ref[i].Id, len(ref[i].Seq))
	}
	words += "@PG\tID:gonomics\tPN:GSW\tVN:1.0\n"
	header.Text = append(header.Text, words)
	return &header
}

func NodesHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		//words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		words = "@SQ\tSN:" + ref[i].Name + "_" + fmt.Sprint(ref[i].Id) + "\tLN:" + fmt.Sprint(len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})
	}
	return &header
}
func indexGenome(genome []*Node, seedLen int) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	var seqCode uint64
	var chromIdx, pos int
	for chromIdx = 0; chromIdx < len(genome); chromIdx++ {

		for pos = 0; pos < len(genome[chromIdx].Seq)-seedLen+1; pos++ {
			seqCode = dnaToNumber(genome[chromIdx].Seq, pos, pos+seedLen)
			answer[seqCode] = append(answer[seqCode], chromAndPosToNumber(chromIdx, pos))
		}
	}
	return answer
}

func MkDictionary(genome []*fasta.Fasta, seedLen int) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	for chromIdx := 0; chromIdx < len(genome); chromIdx++ {
		for pos := 0; pos < len(genome[chromIdx].Seq)-seedLen+1; pos++ {

			seqCode := dnaToNumber(genome[chromIdx].Seq, pos, pos+seedLen)
			answer[seqCode] = append(answer[seqCode], chromAndPosToNumber(chromIdx, pos))
		}
	}
	return answer
}

func chromAndPosToNumber(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func dnaToNumber(seq []dna.Base, start int, end int) uint64 {
	var answer uint64 = uint64(seq[start])
	for i := start + 1; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(seq[i])
	}
	return answer
}

func maskSeedLen(seq []dna.Base, start int, end int, numMismatch int) []uint64 {
	var answer []uint64
	var maskOne uint64 = 1
	var maskTwo uint64 = 2
	var maskThree uint64 = 3
	for i := 1; i < numMismatch; i++ {
		maskOne = maskOne << 1
		maskTwo = maskTwo << 2
		maskThree = maskThree << 3
	}
	answer = append(answer, dnaToNumber(seq, start, end)^maskOne)
	answer = append(answer, dnaToNumber(seq, start, end)^maskTwo)
	answer = append(answer, dnaToNumber(seq, start, end)^maskThree)
	return answer
}

func numberToChromAndPos(code uint64) (int64, int64) {
	var rightSideOnes uint64 = 4294967295
	var leftSideOnes uint64 = rightSideOnes << 32
	var chromIdx uint64 = code & leftSideOnes
	chromIdx = chromIdx >> 32
	var pos uint64 = code & rightSideOnes
	return int64(chromIdx), int64(pos)
}

/*
func SClipCigar(front int64, back int64, lengthOfRead int64, cig []*cigar.Cigar) []*cigar.Cigar {
	var answer []*cigar.Cigar
	if front > 0 {
		answer = append(answer, &cigar.Cigar{RunLength: front, Op: 'S'})
		answer = append(answer, cig...)
		if back < lengthOfRead {
			answer = append(answer, &cigar.Cigar{RunLength: lengthOfRead - back, Op: 'S'})
		}
	} else if back < lengthOfRead {
		answer = append(answer, cig...)
		answer = append(answer, &cigar.Cigar{RunLength: lengthOfRead - back, Op: 'S'})
	} else {
		return cig
	}
	return answer
}*/

func ReadDictionary(filename string) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	file, _ := os.Open(filename)
	defer file.Close()
	reader := bufio.NewReader(file)
	var err error
	var line string
	var words []byte
	for ; err != io.EOF; words, _, err = reader.ReadLine() {
		line = string(words[:])
		data := strings.Split(line, " ")
		for i := 1; i < len(data); i++ {
			answer[common.StringToUint64(data[0])] = append(answer[common.StringToUint64(data[0])], common.StringToUint64(data[i]))
		}
	}
	return answer
}

func WriteDictToFileHandle(file *os.File, input map[uint64][]uint64) error {
	var err error

	for i := range input {
		_, err = fmt.Fprintf(file, "%v", i)
		for j := range input[i] {
			_, err = fmt.Fprintf(file, "\t%v", input[i][j])
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
	}
	return err
}

func WriteDictionary(filename string, data map[uint64][]uint64) {
	//func Write(filename string, data map[int64][]ChrDict) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	WriteDictToFileHandle(file, data)
	//WriteChromDictToFileHandle(file, data)
}

/*
func MapSingleFastq(ref []*Node, chromPosHash map[uint64][]uint64, read *fastq.Fastq, seedLen int, m [][]int64, trace [][]rune) *sam.SamAln {
	//var seedBeds []*bed.Bed
	var seeds []Seed = make([]Seed, 256)
	var codedPositions []uint64
	var bStart, bEnd int64
	var maxScore int64
	var extension int64
	var chrom, pos int64
	var bases int
	var i int
	var subRead, hits int
	var score, bestScore int64 = 0, 0
	var lowRef, lowQuery, highQuery int64
	var alignment []*cigar.Cigar
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: []dna.Base{}, Qual: "", Extra: ""}
	for bases = 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	extension = int64(maxScore/600) + int64(len(read.Seq))
	for subRead = 0; subRead < len(read.Seq)-seedLen+1; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seedLen)]
		for hits = 0; hits < len(codedPositions); hits++ {
			chrom, pos = numberToChromAndPos(codedPositions[hits])
			bStart = pos - extension
			if bStart < 0 {
				bStart = 0
			}
			bEnd = pos + extension
			if bEnd > int64(len(ref[chrom].Seq)) {
				bEnd = int64(len(ref[chrom].Seq))
			}
			//seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
			seeds = addSeed(seeds, chrom, bStart, bEnd)
		}
	}
	//seedBeds = bed.MergeBeds(seedBeds)
	//log.Printf("Length of seed beds, %d", len(seedBeds))
	for i = 0; i < len(seeds); i++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = fmt.Sprintf("%d", seeds[i].Id)
			currBest.Pos = lowRef + seeds[i].Start
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), alignment)
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	seeds = seeds[0:0]
	//seedBeds = [0:0]
	fastq.ReverseComplement(read)
	//extension = int64(maxScore/600) + int64(len(reverse.Seq))
	for subRead = 0; subRead < len(read.Seq)-seedLen+1; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seedLen)]
		for hits = 0; hits < len(codedPositions); hits++ {
			chrom, pos = numberToChromAndPos(codedPositions[hits])
			if chrom >= int64(len(ref)) {
				log.Printf("hashing chrom had index out of bounds b/c of Ns")
			} else {
				bStart = pos - extension
				if bStart < 0 {
					bStart = 0
				}
				bEnd = pos + extension
				if bEnd > int64(len(ref[chrom].Seq)) {
					bEnd = int64(len(ref[chrom].Seq))
				}
				//seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
				seeds = addSeed(seeds, chrom, bStart, bEnd)
			}
		}
	}
	//seedBeds = bed.MergeBeds(seedBeds)
	//log.Printf("Length of reverse seed beds, %d", len(seedBeds))
	for i = 0; i < len(seeds); i++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = fmt.Sprintf("%d", seeds[i].Id)
			currBest.Pos = lowRef + seeds[i].Start
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), alignment)
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
			//fmt.Printf("alignment:\n%s\n", cigar.LocalView(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, alignment, highRef))
		}
	}
	//seeds = seeds[0:0]
	//log.Println(currBest.RName, "\t", currBest.Pos, "\t", read.Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
	return &currBest
}*/
