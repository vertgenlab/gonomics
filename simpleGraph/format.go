package simpleGraph

import (
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
	"strings"
	"sync"
)

func PairedEndAlignFormat(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.PairedSamAln {
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GswFaFormat(gg, readPair.Fwd, seedHash, seedLen, stepSize, m, trace)
	mappedPair.RevSam = GswFaFormat(gg, readPair.Rev, seedHash, seedLen, stepSize, m, trace)
	mappedPair.FwdSam.Flag += 64
	mappedPair.RevSam.Flag += 128
	return &mappedPair
}

func gswPairEndFormat(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.PairedEnd, outgoingSams chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	for read := range incomingFastqs {
		outgoingSams <- PairedEndAlignFormat(gg, read, seedHash, seedLen, stepSize, m, trace)
	}
	wg.Done()
}

func faSize(fa []*fasta.Fasta, fileName string) string {
	var header string = ""
	//header := fmt.Sprintf("##reference=%s\n", fileName)
	for i := 0; i < len(fa); i++ {
		header += fmt.Sprintf("#%s\t%d\n", fa[i].Name, len(fa[i].Seq))
	}
	return header
}

func writeFaHeader(file io.Writer, header string) {
	var err error
	_, err = fmt.Fprintf(file, header)
	common.ExitIfError(err)
}

//TODO: This is the more elegant solution to generate the header, but could not get it to work with samtools...
func DevHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		//words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", ref[i].Name, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})
	}
	return &header
}

func WriteFormatFa(filename string, sg *SimpleGraph, header string) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	writeFaHeader(file, header)
	WriteToGraphHandle(file, sg, lineLength)
}

func ProcessHeader(er *fileio.EasyReader) map[string]*chromInfo.ChromInfo {
	var line string
	var err error
	var nextBytes []byte
	answer := make(map[string]*chromInfo.ChromInfo)
	var words []string
	var i int64
	for nextBytes, err = er.Peek(1); nextBytes[0] == '#' && err == nil; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		words = strings.Split(line, "\t")
		curr := chromInfo.ChromInfo{Name: words[0][1:], Size: common.StringToInt64(words[1]), Order: i}
		answer[words[0][1:]] = &curr
		//answer = append(answer, &chromInfo.ChromInfo{Name: words[0][1:], Size: common.StringToInt64(words[1]), Order: i})
		i++
	}
	return answer
}

func GraphAlignToFaFormat(align *sam.SamAln) {
	if strings.Compare(align.RName, "*") != 0 {
		words := strings.Split(align.RName, "_")
		if len(words) == 3 {
			align.RName = words[0]
			align.Pos += common.StringToInt64(words[2])
		}
	}
}

func GswFaFormat(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, bestScore int64
	var leftPath, rightPath, bestPath []uint32

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
	SortSeedExtended(seeds)
	var tailSeed *SeedDev
	var seedScore int64

	for i = 0; i < len(seeds) && seedCouldBeBetter(seeds[i], bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		tailSeed = toTail(seeds[i])
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}

		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, extension, currRead.Seq[:seeds[i].QueryStart], m, trace)
		seedScore = BlastSeed(seeds[i], currRead)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension, currRead.Seq[tailSeed.QueryStart+tailSeed.Length:], m, trace)

		currScore = leftScore + seedScore + rightScore

		if currScore > bestScore {
			bestPath = CatPaths(CatPaths(leftPath, getSeedPath(seeds[i])), rightPath)
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget + 1)
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

func GswAlignFaFormat(ref *SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, header *sam.SamHeader) {

	//var seedLen int = kMer
	var stepSize int = seedLen - 1
	var numWorkers int = threads
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	//header := devHeader(ref.Nodes)

	var wgAlign, wgWrite sync.WaitGroup

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEnd, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.PairEndToChan(readOne, readTwo, fastqPipe)

	wgAlign.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswPairEndFormat(ref, seedHash, seedLen, stepSize, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}

//TODO: thinking about the idea of making a slice of graphs instead of one graph
/*
func PrintTestGraphFormat(gg []*SimpleGraph) {
	//Write("/dev/stdout", gg)
	lineLength := 50
	file := fileio.EasyCreate("/dev/stdout")
	defer file.Close()
	NodesToFile(file, gg, lineLength)
	EdgesToFile(file, gg)
}

func NodesToFile(file io.Writer, genome []*SimpleGraph, lineLength int) {
	var err error
	var i, j, k int
	for i = 0; i < len(genome); i++ {
		for j = 0; j < len(genome[i].Nodes); j++ {
			_, err = fmt.Fprintf(file, "%s\n", ">"+genome[i].Nodes[j].Name)
			for k = 0; k < len(genome[i].Nodes[j].Seq); k += lineLength {
				if k+lineLength > len(genome[i].Nodes[j].Seq) {
					_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(genome[i].Nodes[j].Seq[k:]))
					common.ExitIfError(err)
				} else {
					_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(genome[i].Nodes[j].Seq[k:k+lineLength]))
					common.ExitIfError(err)
				}
			}
		}
	}
}

func EdgesToFile(file io.Writer, genome []*SimpleGraph) {
	var err error
	var i, j, k int
	for i = 0; i < len(genome); i++ {
		for j = 0; j < len(genome[i].Nodes); j++ {
			_, err = fmt.Fprintf(file, "%s:%d", genome[i].Nodes[j].Name, genome[i].Nodes[j].Id)
			near := genome[i].Nodes[j].Next
			for k = 0; k < len(near); k++ {
				_, err = fmt.Fprintf(file, "\t%v:%v", near[k].Dest.Id, near[k].Prob)
				common.ExitIfError(err)
			}
		_, err = fmt.Fprintf(file, "\n")
		common.ExitIfError(err)
		}
	}
}*/
