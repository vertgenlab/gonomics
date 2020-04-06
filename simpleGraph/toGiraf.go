package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math"
	"path/filepath"
	"strings"
	"sync"
	"time"
)

func GraphSmithWatermanToGiraf(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *giraf.Giraf {
	var currBest giraf.Giraf = giraf.Giraf{
		QName:     read.Name,
		QStart:    0,
		QEnd:      0,
		PosStrand: true,
		Path:      &giraf.Path{},
		Aln:       []*cigar.Cigar{&cigar.Cigar{Op: '*'}},
		AlnScore:  0,
		MapQ:      255,
		Seq:       read.Seq,
		Qual:      []uint8{},
		Notes:     []giraf.Note{giraf.Note{Tag: "XO", Type: 'Z', Value: "~"}},
	}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var minTarget, maxTarget int
	var minQuery, maxQuery int
	var leftScore, rightScore int64 = 0, 0
	var leftPath, rightPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds *SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix, memoryPool)
	SortSeedDevListByTotalLen(&seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev
	for currSeed = seeds; currSeed != nil && seedCouldBeBetterScores(int64(currSeed.TotalLength), int64(currBest.AlnScore), perfectScore, int64(len(read.Seq)), scoreMatrix); currSeed = currSeed.Next {
		tailSeed = getLastPart(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		seedScore = scoreSeedSeq(currSeq, currSeed.QueryStart, tailSeed.QueryStart+tailSeed.Length, scoreMatrix)
		if int(currSeed.TotalLength) == len(currSeq) {
			currScore = seedScore
			minTarget = int(currSeed.TargetStart)
			maxTarget = int(tailSeed.TargetStart + tailSeed.Length)
			minQuery = int(currSeed.QueryStart)
			maxQuery = int(currSeed.TotalLength - 1)
		} else {
			leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[currSeed.TargetId], []dna.Base{}, int(currSeed.TargetStart), []uint32{}, extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], m, trace)
			rightAlignment, rightScore, maxTarget, maxQuery, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension-int(currSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
		}
		currScore = leftScore + seedScore + rightScore
		if currScore > int64(currBest.AlnScore) {
			currBest.QStart = minQuery
			currBest.QEnd = maxQuery
			currBest.PosStrand = currSeed.PosStrand
			currBest.Path = setPath(currBest.Path, minTarget, CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath), maxTarget)
			currBest.Aln = AddSClip(minQuery, len(currSeq), cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(currSeed)), Op: 'M'}), rightAlignment))
			currBest.AlnScore = int(currScore)
			currBest.Seq = currSeq
			if gg.Nodes[currBest.Path.Nodes[0]].Info != nil {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, gg.Nodes[currBest.Path.Nodes[0]].Info.Start)
				currBest.Notes = append(currBest.Notes, infoToNotes(gg.Nodes, currBest.Path.Nodes))
			} else {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, 1)
			}
		}
	}
	if seeds != nil {
		tailSeed = toTail(seeds)
		tailSeed.Next = *memoryPool
		*memoryPool = seeds
	}
	return &currBest
}

func setPath(p *giraf.Path, targetStart int, nodes []uint32, targetEnd int) *giraf.Path {
	p.TStart = targetStart
	p.Nodes = nodes
	p.TEnd = targetEnd
	return p
}

func vInfoToValue(n *Node) string {
	var answer string
	switch {
	case n.Info.Variant == 1:
		answer = fmt.Sprintf("%d=%s", n.Id, "snp")
	case n.Info.Variant == 2:
		answer = fmt.Sprintf("%d=%s", n.Id, "ins")
	case n.Info.Variant == 3:
		answer = fmt.Sprintf("%d=%s", n.Id, "del")
	}
	return answer
}

func infoToNotes(nodes []*Node, path []uint32) giraf.Note {
	var vInfo giraf.Note = giraf.Note{Tag: "XV", Type: 'Z'}
	vInfo.Value = fmt.Sprintf("%d_%d", nodes[0].Info.Allele, nodes[0].Info.Variant)
	if len(path) > 0 {
		for i := 1; i < len(path); i++ {
			if nodes[i].Info.Variant > 0 {
				vInfo.Value += fmt.Sprintf(",%s", vInfoToValue(nodes[path[i]]))
			} else {
				vInfo.Value += fmt.Sprintf(",%d_%d", nodes[i].Info.Allele, nodes[path[i]].Info.Variant)
			}

		}
	}
	return vInfo
}

func routineFqToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.FastqBig, outputChan chan<- *giraf.Giraf, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])
	for read := range inputChan {
		outputChan <- GraphSmithWatermanToGiraf(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, &memStart)
	}
	wg.Done()
}

func GswToGiraf(ref *SimpleGraph, readOne string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64) {
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup

	fastqPipe := make(chan *fastq.FastqBig, 824)
	girafPipe := make(chan *giraf.Giraf, 824)
	go fastq.ReadBigToChan(readOne, fastqPipe)
	log.Printf("Scoring matrix used:\n%s\n", viewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning %s to genome graph...", strings.Split(filepath.Base(readOne), ".")[0])
	start := time.Now()
	for i := 0; i < threads; i++ {
		go routineFqToGiraf(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, girafPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go giraf.GirafChanToFile(output, girafPipe, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(girafPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func WrapPairGiraf(gg *SimpleGraph, readPair *fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{Fwd: nil, Rev: nil}
	mappedPair.Fwd = GraphSmithWatermanToGiraf(gg, readPair.Fwd, seedHash, seedLen, stepSize, scoreMatrix, m, trace, memoryPool)
	mappedPair.Rev = GraphSmithWatermanToGiraf(gg, readPair.Rev, seedHash, seedLen, stepSize, scoreMatrix, m, trace, memoryPool)
	return &mappedPair
}

func routineFqPairToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan *fastq.PairedEndBig, output chan<- *giraf.GirafPair, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])
	for read := range input {
		output <- WrapPairGiraf(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, &memStart)
	}
	wg.Done()
}

func GswToGirafPair(ref *SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64) {
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup

	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	girafPipe := make(chan *giraf.GirafPair, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)
	log.Printf("Scoring matrix used:\n%s\n", viewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)

	log.Printf("Aligning %s and %s to genome graph...", strings.Split(filepath.Base(readOne), ".")[0], strings.Split(filepath.Base(readTwo), ".")[0])
	start := time.Now()
	for i := 0; i < threads; i++ {
		go routineFqPairToGiraf(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, girafPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go giraf.GirafPairChanToFile(output, girafPipe, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(girafPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func WrapGirafLiftoverToSam(ref *SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64, header *sam.SamHeader) {
	log.SetFlags(log.Ldate | log.Ltime)
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	//log.Printf("Setting up read and write channels...\n\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)

	log.Printf("Scoring matrix used:\n%s\n", viewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning sequence to genome graph...")
	start := time.Now()
	for i := 0; i < threads; i++ {
		go routineGirafToSam(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func GirafToSam(ag *giraf.Giraf, read *fastq.FastqBig) *sam.SamAln {
	curr := &sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0\tGP:Z:-1\tXO:Z:~"}
	if strings.Compare(ag.Notes[0].Value, "~") == 0 {
		//read is unMapped
		return curr
	} else {
		target := strings.Split(ag.Notes[0].Value, "=")
		curr.RName = target[0]
		curr.Pos = int64(ag.Path.TStart) + common.StringToInt64(target[1])
		curr.Flag = getSamFlags(ag)
		curr.Cigar = ag.Aln
		curr.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:Z:%d\t%s", ag.AlnScore, PathToString(ag.Path.Nodes), ag.Path.TStart, giraf.NoteToString(ag.Notes[1]))
	}
	return curr
}

func GirafLiftoverToSam(gg *SimpleGraph, readPair *fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *sam.PairedSamAln {
	mappedPair := WrapPairGiraf(gg, readPair, seedHash, seedLen, stepSize, scoreMatrix, m, trace, memoryPool)
	toSamPair := GirafPairToSam(mappedPair, readPair)
	return toSamPair
}

func GirafPairToSam(ag *giraf.GirafPair, readPair *fastq.PairedEndBig) *sam.PairedSamAln {
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: &sam.SamAln{}, RevSam: &sam.SamAln{}}
	mappedPair.FwdSam = GirafToSam(ag.Fwd, readPair.Fwd)
	mappedPair.RevSam = GirafToSam(ag.Rev, readPair.Rev)
	mappedPair.FwdSam.Flag += 64
	mappedPair.RevSam.Flag += 128
	if isProperPairAlign(ag) {
		mappedPair.FwdSam.Flag += 2
		mappedPair.RevSam.Flag += 2
	}
	return &mappedPair
}

func routineGirafToSam(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.PairedEndBig, outputChan chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])
	for read := range inputChan {
		outputChan <- GirafLiftoverToSam(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, &memStart)
	}
	wg.Done()
}

func isProperPairAlign(mappedPair *giraf.GirafPair) bool {
	if math.Abs(float64(mappedPair.Fwd.Path.TStart-mappedPair.Rev.Path.TStart)) < 10000 {
		if mappedPair.Fwd.Path.TStart < mappedPair.Rev.Path.TStart && mappedPair.Fwd.PosStrand && !mappedPair.Rev.PosStrand {
			return true
		}
		if mappedPair.Fwd.Path.TStart > mappedPair.Rev.Path.TStart && !mappedPair.Fwd.PosStrand && mappedPair.Rev.PosStrand {
			return true
		}
	}
	return false
}

func getSamFlags(ag *giraf.Giraf) int64 {
	var answer int64
	if !ag.PosStrand {
		answer += 16
	}
	if ag.AlnScore < 1200 {
		answer += 4
	}
	return answer
}
