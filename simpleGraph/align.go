package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/chromInfo"
	"strconv"
	"runtime"
	"sync"
	"log"
	"os"
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

func chromAndPosToNumber(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func dnaToNumber(dna []dna.Base, start int, end int) uint64 {
	var answer uint64 = uint64(dna[start])
	for i := start + 1; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(dna[i])
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

func MapFastq(ref []*Node, read *fastq.Fastq, seed int, chromPosHash map[uint64][]uint64, c chan *sam.SamAln) {
	//chromPosHash := indexGenome(ref, seed)
	//var answer []*sam.SamAln
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
	//seedBeds = nil
	for bases = 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	//fmt.Printf("Max Score: %d\n", maxScore)
	extension = int64(maxScore/600) + int64(len(read.Seq))
	for subRead = 0; subRead < len(read.Seq)-seed+1; subRead++ {

		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seed)]
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
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], read.Seq, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = seedBeds[beds].Chrom
			//currBest.RName = bestBed.Chrom
			currBest.Pos = lowRef + seedBeds[beds].ChromStart
			//currBest.Pos = lowRef+bestBed.ChromStart
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), cigar.FromString(align.PrintCigar(alignment)))
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	//reverse read
	seedBeds = nil
	//dna.ReverseComplement(read.Seq)
	reverse := fastq.ReverseComplementFastq(read)
	//for bases = 0; bases < len(reverse.Seq); bases++ {
	//	maxScore += HumanChimpTwoScoreMatrix[reverse.Seq[bases]][reverse.Seq[bases]]
	//}
	//fmt.Printf("Max Score: %d\n", maxScore)
	extension = int64(maxScore/600) + int64(len(reverse.Seq))
	for subRead = 0; subRead < len(reverse.Seq)-seed+1; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(reverse.Seq, subRead, subRead+seed)]
		for hits = 0; hits < len(codedPositions); hits++ {
			chrom, pos = numberToChromAndPos(codedPositions[hits])
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
		}
	}
	seedBeds = bed.MergeBeds(seedBeds)
	for beds = 0; beds < len(seedBeds); beds++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], reverse.Seq, HumanChimpTwoScoreMatrix, -600)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = seedBeds[beds].Chrom
			//currBest.RName = bestBed.Chrom
			currBest.Pos = lowRef + seedBeds[beds].ChromStart
			//currBest.Pos = lowRef+bestBed.ChromStart
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(reverse.Seq)), cigar.FromString(align.PrintCigar(alignment)))
			currBest.Seq = reverse.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
		currBest.Cigar = []*cigar.Cigar{&cigar.Cigar{RunLength: 0, Op: '*'}}
	}
	log.Println(currBest.RName, "\t", read.Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
	c <- &currBest
}

func GSW(ref []*Node, m map[uint64][]uint64, fastqFile string, samFile string) {
	var answer *sam.SamAln
	header := QFragHeader(ref)
	
	outFile, _ := os.Create(samFile+".sam")
	defer outFile.Close()

	sam.WriteHeaderToFileHandle(outFile, header)

	inputCh := make(chan *fastq.Fastq)
	outputCh := make(chan *sam.SamAln)
	
	var wg sync.WaitGroup
	runtime.GOMAXPROCS(runtime.NumCPU())
	threads := runtime.NumCPU()*50
	wg.Add(threads)
	fmt.Println("Num of threads: ", threads)
	//var wg sync.WaitGroup
	//wg.Add(100)
	var j int
	file := fileio.EasyOpen(fastqFile)
	defer file.Close()
	var done bool
	
	var fq *fastq.Fastq
	var countReads int = 0
	var ok bool
	var input *fastq.Fastq
	
	for workers := 0; workers < threads; workers++ {
		go func() {
			for {
				input, ok = <-inputCh
				if !ok {
					wg.Done()
					return
				} else {
					MapFastq(ref, input, 25, m, outputCh)
				}
			}
		}()
	}

	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		//go gsw(ref, reads[i], c)
		countReads++
		inputCh <-fq
		//wg.Add(1)
		//go warrior(ref, fq, 20, m, c)
	}
	close(inputCh)
	for j = 0; j < countReads; j++ {
		answer = <-outputCh
		sam.WriteAlnToFileHandle(outFile, answer)
	}

	wg.Wait()
}

func MapReads(ref []*Node, reads []*fasta.Fasta, seed int) []*sam.SamAln {
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

	var score, bestScore int64 = 0, 0
	var lowRef, lowQuery, highQuery int64
	//initial sam alignment pointer
	var alignment []align.Cigar
	for i := 0; i < len(reads); i++ {
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
		//if bestScore < 1200 {
		//	currBest.Flag = 4
		//	currBest.Cigar = []*cigar.Cigar{&cigar.Cigar{RunLength: 0, Op: '*'}}
		//}
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

func BasicHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted") 
	var words string

	for i := 0; i < len(ref); i++ {
		words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: strconv.FormatInt(ref[i].Id, 10), Size: int64(len(ref[i].Seq))})
		
	}
	return &header
}
