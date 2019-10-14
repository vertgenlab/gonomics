package simpleGraph

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
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
	"runtime"
	"strconv"
	"strings"
	"sync"
)

func indexGenome(genome []*Node, seedLen int) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	for chromIdx := 0; chromIdx < len(genome); chromIdx++ {
		for pos := 0; pos < len(genome[chromIdx].Seq)-seedLen+1; pos++ {
			seqCode := dnaToNumber(genome[chromIdx].Seq, pos, pos+seedLen)
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

func MapFastq(ref []*fasta.Fasta, read *fastq.Fastq, seed int, chromPosHash map[uint64][]uint64, file *os.File) {
//func MapFastq(ref []*fasta.Fasta, read *fastq.Fastq, seed int, chromPosHash map[uint64][]uint64, file *os.File) {
	//defer wg.Done()
	//chromPosHash := indexGenome(ref, seed)
	//var answer []*sam.SamAln
	fasta.AllToUpper(ref)
	var seedBeds []*bed.Bed
	var codedPositions []uint64
	var bStart, bEnd int64
	var maxScore int64
	var extension int64
	var chrom, pos int64
	//calculate max score
	var bases, beds int
	var subRead, hits int

	var score, bestScore int64 = 0, 0
	var lowRef, lowQuery, highQuery int64
	//initial sam alignment pointer
	var alignment []align.Cigar

	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: []dna.Base{}, Qual: "", Extra: ""}

	for bases = 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}

	extension = int64(maxScore/600) + int64(len(read.Seq))
	for subRead = 0; subRead < len(read.Seq)-seed; subRead++ {

		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seed)]
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
				//if bEnd > int64(len(chrSize[m[putTogether(read.Seq[i:i+seed])][j].Chr])) {
				if bEnd > int64(len(ref[chrom].Seq)) {
					bEnd = int64(len(ref[chrom].Seq))
				}
				seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
			}
		}

	}
	seedBeds = bed.MergeBeds(seedBeds)
	//seedBeds = bed.HighestScoreBed(seedBeds)
	for beds = 0; beds < len(seedBeds); beds++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], read.Seq, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = ref[common.StringToInt64(seedBeds[beds].Chrom)].Name
			//currBest.RName = bestBed.Chrom
			currBest.Pos = lowRef + seedBeds[beds].ChromStart
			//currBest.Pos = lowRef+bestBed.ChromStart
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), cigar.FromString(align.PrintCigar(alignment)))
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	seedBeds = nil
	reverse := fastq.ReverseComplementFastq(read)
	extension = int64(maxScore/600) + int64(len(reverse.Seq))
	for subRead = 0; subRead < len(reverse.Seq)-seed; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(reverse.Seq, subRead, subRead+seed)]
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
				seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
			}

		}
	}
	seedBeds = bed.MergeBeds(seedBeds)
	//seedBeds = bed.HighestScoreBed(seedBeds)
	for beds = 0; beds < len(seedBeds); beds++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], reverse.Seq, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = ref[common.StringToInt64(seedBeds[beds].Chrom)].Name
			//currBest.RName = bestBed.Chrom
			currBest.Pos = lowRef + seedBeds[beds].ChromStart
			//currBest.Pos = lowRef+bestBed.ChromStart
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(reverse.Seq)), cigar.FromString(align.PrintCigar(alignment)))
			currBest.Seq = reverse.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	//if bestScore < 1200 {
	//	currBest.Flag = 4
	//	currBest.Cigar = []*cigar.Cigar{&cigar.Cigar{RunLength: 0, Op: '*'}}
	//}
	log.Println(currBest.RName, "\t", currBest.Pos, "\t", read.Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
	//sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: []dna.Base{}, Qual: "", Extra: ""}
	sam.WriteAlnToFileHandle(file, &currBest)
	//c <- &currBest
	//wg.Done()
}
/*
func GSW(ref []*fasta.Fasta, m map[uint64][]uint64, fastqFile string, samFile string) {
	//var answer *sam.SamAln
	header := BasicHeader(ref)
	wg := new(sync.WaitGroup)
	outFile, _ := os.Create(samFile)
	defer outFile.Close()

	sam.WriteHeaderToFileHandle(outFile, header)

	//c := make(chan *sam.SamAln)

	file := fileio.EasyOpen(fastqFile)
	defer file.Close()
	var done bool

	var fq *fastq.Fastq
	//var countReads int = 0

	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		wg.Add(1)
		go MapFastq(ref, fq, 25, m, wg, outFile)
	}
	wg.Wait()
}*/

func MapReads(ref []*Node, reads []*fastq.Fastq, seed int) []*sam.SamAln {
	chromPosHash := indexGenome(ref, seed)
	var answer []*sam.SamAln
	var seedBeds []*bed.Bed
	var codedPositions []uint64
	var bStart, bEnd int64
	var maxScore int64
	var extension int64
	var chrom, pos int64
	//calculate max score
	var bases, beds int
	var subRead, hits int

	var score, bestScore int64
	var lowRef, lowQuery, highQuery int64
	//initial sam alignment pointer
	var alignment []align.Cigar
	for i := 0; i < len(reads); i++ {
		score, bestScore = 0, 0
		var currBest sam.SamAln = sam.SamAln{QName: reads[i].Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: reads[i].Seq, Qual: "", Extra: ""}
		seedBeds = nil
		for bases = 0; bases < len(reads[i].Seq); bases++ {
			maxScore += HumanChimpTwoScoreMatrix[reads[i].Seq[bases]][reads[i].Seq[bases]]
		}
		//fmt.Printf("Max Score: %d\n", maxScore)
		extension = int64(maxScore/600) + int64(len(reads[i].Seq))

		for subRead = 0; subRead < len(reads[i].Seq)-seed+1; subRead++ {

			codedPositions = chromPosHash[dnaToNumber(reads[i].Seq, subRead, subRead+seed)]
			for hits = 0; hits < len(codedPositions); hits++ {
				chrom, pos = numberToChromAndPos(codedPositions[hits])
				//fmt.Printf("Chrom: %d Position: %d\n", chrom, pos)
				//fmt.Printf("Length of hits: %d\n", len(codedPositions))
				bStart = pos - extension
				if bStart < 0 {
					bStart = 0
				}
				bEnd = pos + extension
				//if bEnd > int64(len(chrSize[m[putTogether(read.Seq[i:i+seed])][j].Chr])) {
				if bEnd > int64(len(ref[chrom].Seq)) {
					bEnd = int64(len(ref[chrom].Seq))
				}
				seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(ref[chrom].Id, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
				//fmt.Println(bed.Bed{Chrom: strconv.FormatInt(ref[chrom].Id, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
			}

		}
		seedBeds = bed.MergeBeds(seedBeds)
		for beds = 0; beds < len(seedBeds); beds++ {
			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], reads[i].Seq, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score
				currBest.Flag = 0
				currBest.RName = seedBeds[beds].Chrom
				//currBest.RName = bestBed.Chrom
				currBest.Pos = lowRef + seedBeds[beds].ChromStart
				//currBest.Pos = lowRef+bestBed.ChromStart
				currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(reads[i].Seq)), cigar.FromString(align.PrintCigar(alignment)))
				currBest.Seq = reads[i].Seq
				//currBest.Qual = string(read.Qual)

			}
		}
		if bestScore < 1200 {
			currBest.Flag = 4
			currBest.Cigar = []*cigar.Cigar{&cigar.Cigar{RunLength: 0, Op: '*'}}
		}
		//fmt.Println(currBest)
		fmt.Println(currBest.RName, "\t", reads[i].Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
		answer = append(answer, &currBest)
	}

	return answer
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

func SClipCigar(front int64, back int64, lengthOfRead int64, cig []*cigar.Cigar) []*cigar.Cigar {
	var answer []*cigar.Cigar
	if front > 0 {
		answer = append(answer, &cigar.Cigar{RunLength: front, Op: 'S'})
		answer = append(answer, cig...)
		if back < lengthOfRead {
			//add soft clip on both back and front
			//Cigar = &Cigar{RunLength: front, Op: 'S'}

			answer = append(answer, &cigar.Cigar{RunLength: lengthOfRead - back, Op: 'S'})
		}
	} else if back < lengthOfRead {
		answer = append(answer, cig...)
		answer = append(answer, &cigar.Cigar{RunLength: lengthOfRead - back, Op: 'S'})
	} else {
		answer = cig
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

func GSW(ref []*fasta.Fasta, m map[uint64][]uint64, fastqFile string, samFile string) {
	//var query *fastq.Fastq
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

	for workers := 0; workers < threads; workers++ {
		go func(jobNumber int) {
			defer wg.Done()
			for {
				query, ok := <-tasks
				if !ok {
					return
				}
				MapFastq(ref, query, 25, m, outFile)
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
	//MapFastq(ref, j, 25, m, wg, outFile)
	
	//close(c)
	//wg.Wait()
}
