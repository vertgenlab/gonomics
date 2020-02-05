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
	"log"
	"os"
	"strconv"
	"strings"
)

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
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

func numberToChromAndPos(code uint64) (int64, int64) {
	var rightSideOnes uint64 = 4294967295
	var leftSideOnes uint64 = rightSideOnes << 32
	var chromIdx uint64 = code & leftSideOnes
	chromIdx = chromIdx >> 32
	var pos uint64 = code & rightSideOnes
	return int64(chromIdx), int64(pos)
}

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
}

func scoreSeed(seed *SeedDev, read *fastq.Fastq) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

//func GraphSmithWaterman(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
func GraphSmithWaterman(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, seedScore, bestScore int64
	var leftPath, rightPath, bestPath []uint32

	var currScore, maxScore int64 = 0, 0
	for bases := 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	ext := int(maxScore/600) + len(read.Seq)
	m, trace := SwMatrixSetup(int64(ext + 1))

	var currRead *fastq.Fastq = nil
	var seeds []*SeedDev = findSeedsInMap(seedHash, read, seedLen, stepSize, true)
	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInMap(seedHash, revCompRead, seedLen, stepSize, false)
	seeds = append(seeds, revCompSeeds...)
	SortSeedDevByLen(seeds)
	//log.Printf("number of seed hits = %d\n", len(seeds))
	if len(seeds) > 24 {
		seeds = seeds[:8]
		//log.Printf("number of seed hits changed: %d\n", len(seeds))
	}
	for i = 0; i < len(seeds) && seedCouldBeBetter(seeds[i], bestScore, maxScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}
		seedScore = scoreSeed(seeds[i], currRead)
		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, ext, currRead.Seq[:seeds[i].QueryStart], m, trace)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart+seeds[i].Length), []uint32{}, ext, currRead.Seq[seeds[i].QueryStart+seeds[i].Length:], m, trace)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 1
			}
			bestPath = CatPaths(AddPath(seeds[i].TargetId, leftPath), rightPath)
			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget)
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(seeds[i].Length), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(int64(minQuery), int64(len(currRead.Seq)), currBest.Cigar)
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(bestPath, gg)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
		currBest.Cigar = []*cigar.Cigar{&cigar.Cigar{RunLength: 0, Op: '*'}}
	}
	return &currBest
}

func AddSClip(front int64, lengthOfRead int64, cig []*cigar.Cigar) []*cigar.Cigar {
	back := front + cigar.QueryLength(cig)
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
		//answer = cig
		return cig
	}
	return answer
}

func MaxScore(seq []dna.Base, scores [][]int64) int64 {
	var maxScore int64 = 0
	for i := 0; i < len(seq); i++ {
		maxScore += scores[seq[i]][seq[i]]
	}
	return maxScore
}

/*
func GSW(ref []*fasta.Fasta, m map[uint64][]uint64, seedLen int, fastqFile string, samFile string) {
	//var answer *sam.SamAln
	header := BasicHeader(ref)
	wg := new(sync.WaitGroup)
	outFile, _ := os.Create(samFile)
	defer outFile.Close()

	sam.WriteHeaderToFileHandle(outFile, header)
	runtime.GOMAXPROCS(runtime.NumCPU())
	//c := make(chan *sam.SamAln)

	file := fileio.EasyOpen(fastqFile)
	defer file.Close()
	var done bool

	var fq *fastq.Fastq
	//var countReads int = 0
	mSet, trace := swMatrixSetup(10000)
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		wg.Add(1)
		go MapFastq(ref, fq, seedLen, m, mSet, trace, outFile)
	}
	wg.Wait()
}*/

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
		//answer = cig
		return cig
	}
	return answer
}

func BasicHeader(ref []*fasta.Fasta) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		//words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		words = "@SQ\tSN:" + ref[i].Name + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})

	}
	return &header
}

func NodesHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		//words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		words = "@SQ\tSN:" + ref[i].Name + "_" + fmt.Sprint(ref[i].Id) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})

	}
	return &header
}

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
		//fmt.Println(data)
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

//call aligner with goroutines
/*
func GSW(ref []*fasta.Fasta, m map[uint64][]uint64, seedLen int, fastqFile string, samFile string) {
	fasta.AllToUpper(ref)
	header := BasicHeader(ref)
	wg := new(sync.WaitGroup)
	outFile, _ := os.Create(samFile)
	defer outFile.Close()
	sam.WriteHeaderToFileHandle(outFile, header)
	file := fileio.EasyOpen(fastqFile)
	defer file.Close()

	runtime.GOMAXPROCS(runtime.NumCPU())
	//threads := runtime.NumCPU()
	threads := runtime.NumCPU()
	tasks := make(chan *fastq.Fastq)
	mSet, trace := swMatrixSetup(10000)
	for workers := 0; workers < threads; workers++ {
		go func(jobNumber int) {
			defer wg.Done()
			for {
				query, ok := <-tasks
				if !ok {
					return
				}
				MapFastq(ref, query, seedLen, m, mSet, trace, outFile)
			}

		}(workers)
	}

	var fq *fastq.Fastq
	var done bool
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		wg.Add(1)
		tasks <- fq
	}
	close(tasks)
	wg.Wait()
}*/

/*
func GSW(ref []*fasta.Fasta, m map[uint64][]uint64, fastqFile string, samFile string) {
	fasta.AllToUpper(ref)
	header := BasicHeader(ref)

	outFile, _ := os.Create(samFile)
	defer outFile.Close()
	sam.WriteHeaderToFileHandle(outFile, header)
	file := fileio.EasyOpen(fastqFile)
	defer file.Close()

	var fq *fastq.Fastq
	mSet, trace := swMatrixSetup(10000)
	var done bool
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		MapFastq(ref, fq, 25, m, mSet, trace, outFile)
	}

}*/
