package simpleGraph

import(
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/edotau/simpleio"
	"log"
	"strings"
	"bytes"
	"sync"

)

func SimplyGsw(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]simpleio.CigarOp) *giraf.Giraf {
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
		Qual:      read.Qual,
		Notes:     []giraf.Note{giraf.Note{Tag: "XO", Type: 'Z', Value: "~"}},
	
		ByteCigar:       []simpleio.ByteCigar{},

	}
	var leftAlignment, rightAlignment []simpleio.ByteCigar = []simpleio.ByteCigar{}, []simpleio.ByteCigar{}
	var minTarget, maxTarget int
	var minQuery, maxQuery int
	var leftScore, rightScore int64 = 0, 0
	var leftPath, rightPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds []*SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix)
	SortSeedDevByLen(seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev

	var leftSeq, rightSeq []dna.Base

	for i := 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].TotalLength), int64(currBest.AlnScore), perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		currSeed = seeds[i]
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
			leftSeq = make([]dna.Base, 0, currSeed.QueryStart)
			leftAlignment, leftScore, minTarget, minQuery, leftPath = LeftAlignTraversal(gg.Nodes[currSeed.TargetId], leftSeq, int(currSeed.TargetStart), leftPath, extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], m, trace)
			

			rightSeq = make([]dna.Base, 0, uint32(len(currSeq))-tailSeed.QueryStart-tailSeed.Length)
			rightAlignment, rightScore, maxTarget, maxQuery, rightPath = RightAlignTraversal(gg.Nodes[tailSeed.TargetId], rightSeq, int(tailSeed.TargetStart+tailSeed.Length), rightPath, extension-int(currSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
			//log.Printf("left alignment: %s\n", simpleio.ByteCigarString(leftAlignment))
			simpleio.ByteCigarString(rightAlignment)
			//log.Printf("right alignment: %s\n", simpleio.ByteCigarString(rightAlignment))
		}
		currScore = leftScore + seedScore + rightScore
		if currScore > int64(currBest.AlnScore) {
			currBest.QStart = minQuery
			currBest.QEnd = maxQuery
			currBest.PosStrand = currSeed.PosStrand
			
			currBest.Path = setPath(currBest.Path, minTarget, CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath), maxTarget)
			//currBest.Aln = AddSClip(minQuery, len(currSeq), cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(currSeed)), Op: 'M'}), rightAlignment))
			
			currBest.ByteCigar = simpleio.SoftClipBases(minQuery, len(currSeq), simpleio.CatByteCigar(simpleio.AddCigar(leftAlignment, simpleio.ByteCigar{RunLen: sumLen(currSeed), Op: 'M'}), rightAlignment))
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
	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return &currBest
}

func RoutineSimplyGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.PairedEndBig, outputChan chan<- *giraf.GirafPair, wg *sync.WaitGroup) {
	m, trace := simpleio.MatrixSetup(10000)
	for read := range inputChan {
		outputChan <- WrapSimplyGsw(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	}
	wg.Done()
}

func WrapSimplyGsw(gg *SimpleGraph, readPair *fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]simpleio.CigarOp) *giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{Fwd: nil, Rev: nil}
	mappedPair.Fwd = SimplyGsw(gg, readPair.Fwd, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	mappedPair.Rev = SimplyGsw(gg, readPair.Rev, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	setGirafFlags(&mappedPair)
	return &mappedPair
}

func SoftClipBases(front int, lengthOfRead int, cig []simpleio.ByteCigar) []simpleio.ByteCigar {
	var runLen int = simpleio.QueryRunLen(cig)
	if runLen < lengthOfRead {
		answer := make([]simpleio.ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, simpleio.ByteCigar{RunLen: uint32(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+simpleio.QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, simpleio.ByteCigar{RunLen: uint32(lengthOfRead-front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}


}

func SimpleSequenceGraph(ref <-chan *fasta.Fasta, vcfMap map[string][]*vcf.Vcf) *SimpleGraph {
	gg := NewGraph()
	var filterVcf []*vcf.Vcf = make([]*vcf.Vcf, 0)
	for chr := range ref {
		filterVcf = vcfMap[chr.Name]
		if len(filterVcf) != 0 {
			vcf.Sort(filterVcf)
			gg = vChrGraph(gg, chr, filterVcf)
		} else {
			// Given the input is a vcf containing large structural variance
			// It is possible for a chromosome to contain no call variants (ex. chrM).
			// In such a case, we add the entire chromosome as one node.
			chrNode := &Node{Id: uint32(len(gg.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(chr.Seq), Prev: nil, Next: nil, Info: &Annotation{Allele: 0, Start: 1, Variant: 0}}
			AddNode(gg, chrNode)
		}
	}
	return gg
}

func SimpleChrom(genome *SimpleGraph, chr *fasta.Fasta, vcfsChr []*vcf.Vcf) *SimpleGraph {
	vcfsChr = append(vcfsChr, &vcf.Vcf{Chr: chr.Name, Pos: int64(len(chr.Seq))})
	//log.Printf("Found %d variants on %s", len(vcfsChr), chr.Name)
	fasta.ToUpper(chr)
	var currMatch *Node = &Node{}
	var lastMatch *Node = &Node{}
	var refAllele, altAllele *Node = &Node{}, &Node{}
	var prev []*Node = nil
	var i, j, edge int
	var index int64 = 0
	for i = 0; i < len(vcfsChr)-1; i++ {
		if strings.Compare(chr.Name, vcfsChr[i].Chr) != 0 {
			log.Fatalf("Error: chromosome names do not match...\n")
		}
		if vcfsChr[i].Pos-index > 0 {
			currMatch = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(chr.Seq[index : vcfsChr[i].Pos-1]), Prev: nil, Next: make([]*Edge, 0, 2), Info: &Annotation{Allele: 0, Start: uint32(index + 1), Variant: 0}}
			if len(currMatch.Seq) == 0 {
				currMatch = lastMatch
				//assuming we already created the ref allele, we only need to create
				//the alt alleles this iteration
				if vcf.Snp(vcfsChr[i]) {
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)
				} else if vcf.Ins(vcfsChr[i]) {
					insertion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)[1:]), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 2}}
					AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)

					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					deletion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Ref)[1:]), Prev: nil, Next: nil, Info: &Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 3}}
					AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = common.MinInt64(vcfsChr[i].Pos+int64(len(deletion.Seq))-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + int64(len(deletion.Seq))
					}
				} else if strings.Compare(vcfsChr[i].Format, "SVTYPE=SNP;INS") == 0 || strings.Compare(vcfsChr[i].Format, "SVTYPE=SNP;DEL") == 0 || strings.Compare(vcfsChr[i].Format, "SVTYPE=HAP") == 0 {
					altAllele := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = vcfsChr[i].Pos + int64(len(refAllele.Seq)) - 1
				}
				lastMatch = currMatch
			} else if len(currMatch.Seq) > 0 {

				AddNode(genome, currMatch)
				if lastMatch != nil {
					if len(lastMatch.Next) > 0 {
						for edge = 0; edge < len(lastMatch.Next); edge++ {
							AddEdge(lastMatch.Next[edge].Dest, currMatch, 1)
						}
					}
					if i > 0 {
						if vcf.Snp(vcfsChr[i-1]) || isHaplotypeBlock(vcfsChr[i-1]) {
							AddEdge(altAllele, currMatch, 1)
						}
					}
					AddEdge(lastMatch, currMatch, 1)
					SetEvenWeights(lastMatch)
				}
				if vcf.Snp(vcfsChr[i]) {

					refAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Ref)), Prev: nil, Next: nil, Info: &Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 0.5)

					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 1}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 0.5)

					currMatch = refAllele
					index = vcfsChr[i].Pos
					for j = i + 1; j < len(vcfsChr)-1; j++ {
						if vcf.Snp(vcfsChr[j-1]) && vcf.Snp(vcfsChr[j]) && vcfsChr[j].Pos-1 == vcfsChr[j-1].Pos {
							refAllele.Seq = append(refAllele.Seq, dna.StringToBases(vcfsChr[j].Ref)...)
							altAllele.Seq = append(altAllele.Seq, dna.StringToBases(vcfsChr[j].Alt)...)
							index = vcfsChr[j].Pos
						} else {
							lastMatch = currMatch
							i = j - 1
							break
						}
					}
				} else if vcf.Ins(vcfsChr[i]) {
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					insertion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 2}}
					AddNode(genome, insertion)
					AddEdge(currMatch, insertion, 1)
					index = vcfsChr[i].Pos - 1
				} else if vcf.Del(vcfsChr[i]) {
					//TODO: does not deal with ends of large deletions well, might contain extra sequence towards the end.
					//currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Alt)...)
					deletion := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Ref)), Prev: nil, Next: nil, Info: &Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 3}}
					AddNode(genome, deletion)
					AddEdge(currMatch, deletion, 1)
					if strings.Contains(vcfsChr[i].Id, "pbsv") {
						index = common.MinInt64(vcfsChr[i].Pos+int64(len(deletion.Seq))-1, vcfsChr[i+1].Pos-1)
					} else {
						index = vcfsChr[i].Pos + int64(len(deletion.Seq))
					}
				} else if isINV(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					inversion := mkInversionNode(genome, vcfsChr[i], chr)
					AddNode(genome, inversion)
					AddEdge(currMatch, inversion, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isCNV(vcfsChr[i]) || isDup(vcfsChr[i]) {
					currMatch.Seq = append(currMatch.Seq, dna.StringToBases(vcfsChr[i].Ref)...)
					copyVariance := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(chr.Seq[vcfsChr[i].Pos:getSvEnd(vcfsChr[i])]), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 6}}
					AddNode(genome, copyVariance)
					AddEdge(currMatch, copyVariance, 1)
					index = getSvEnd(vcfsChr[i])
				} else if isHaplotypeBlock(vcfsChr[i]) {
					refAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Ref)), Prev: nil, Next: nil, Info: &Annotation{Allele: 0, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, refAllele)
					AddEdge(currMatch, refAllele, 1)
					altAllele = &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(vcfsChr[i].Alt)), Prev: nil, Next: nil, Info: &Annotation{Allele: 1, Start: uint32(vcfsChr[i].Pos), Variant: 4}}
					AddNode(genome, altAllele)
					AddEdge(currMatch, altAllele, 1)
					index = common.MinInt64(vcfsChr[i].Pos+int64(len(refAllele.Seq))-1, vcfsChr[i+1].Pos-1)
					currMatch = refAllele
				}
				lastMatch = currMatch
			} else {
				//num++
				//log.Printf("index=%d, vcfPos=%d\ndiff=%d\n%s", index, vcfsChr[i].Pos, vcfsChr[i].Pos-index, vcfsChr[i].Info)
			}
		}
	}
	//Case: last node
	lastNode := &Node{Id: uint32(len(genome.Nodes)), Name: chr.Name, SeqTwoBit: dnaTwoBit.NewTwoBit(chr.Seq[index:]), Prev: make([]*Edge, 0, len(prev)), Next: nil, Info: &Annotation{Allele: 0, Start: uint32(index + 1), Variant: 0}}
	AddNode(genome, lastNode)
	for edge = 0; edge < len(lastMatch.Next); edge++ {
		AddEdge(lastMatch.Next[edge].Dest, lastNode, 1)
	}
	if vcf.Snp(vcfsChr[len(vcfsChr)-2]) || isHaplotypeBlock(vcfsChr[len(vcfsChr)-2]) {
		AddEdge(altAllele, lastNode, 1)
	}
	AddEdge(lastMatch, lastNode, 1)
	SetEvenWeights(lastMatch)
	return genome
}

func checkAlignment(aln *giraf.Giraf, genome *SimpleGraph) bool {
	qName := strings.Split(aln.QName, "_")
	//if len(qName) < 5 {
	//	log.Fatalf("Error: input giraf file does not match simulation format...\n")
	//}
	if len(aln.Aln) < 1 {
		return false
	}


	
	targetStart := aln.Path.TStart
	targetEnd := aln.Path.TEnd
	//if len(aln.Aln) < 1 {
	if aln.ByteCigar[0].Op == 'S' {
	//log.Printf("%s\n", giraf.GirafToString(aln))
		targetStart = targetStart - int(aln.ByteCigar[0].RunLen)
	}
	if aln.Aln[len(aln.Aln)-1].Op =='S' {
		targetEnd = targetEnd + int(aln.ByteCigar[len(aln.ByteCigar)-1].RunLen)

	//}
	
	}
	if common.StringToInt(qName[0])  == int(aln.Path.Nodes[0]) && common.StringToInt(qName[1]) ==  targetStart && targetEnd == common.StringToInt(qName[3]) {
		//log.Printf("%s\n", giraf.GirafToString(aln))
		//log.Printf("Results: %d != %d or %d != %d\n", headNode, aln.Path.Nodes[0], startPos, aln.Path.TStart)
	//	log.Printf("%s\n", giraf.GirafToString(aln))
		return true
	} else {
		//log.Printf("endPos=%d, right side cigar runLength: %d\n", endPos, aln.Aln[len(aln.Aln)-1].RunLen)
		//log.Printf("%s\n", giraf.GirafToString(aln))
		//log.Printf("Error: this read is not aligning correctly...\n")
	}
return false
}

func isGirafPairCorrect(input <-chan *giraf.GirafPair, genome *SimpleGraph, wg *sync.WaitGroup, numReads int) {
	var unmapped int = 0
	for pair := range input {
		if pair.Fwd.Aln == nil {
			log.Printf("Warning... read was recorded as unmapped...\n")
			continue
		}
		if !checkAlignment(pair.Fwd, genome) {
			unmapped++
			//log.Printf("Error: failed alignment simulation...\n")
			//log.Printf("%s\n", giraf.GirafToString(pair.Fwd))
			//log.Printf("%s\n", giraf.GirafToString(pair.Rev))	
		}
		if !checkAlignment(pair.Rev, genome) {
			//log.Printf("Error: failed alignment simulation...\n")
			//log.Printf("%s\n", giraf.GirafToString(pair.Rev))
			unmapped++
		}
	}
	
	log.Printf("Mapped %d out of %d\n", numReads-unmapped, numReads)
	log.Printf("%f of the reads are mapping correctly\n", percentOfFloat(numReads-unmapped, numReads))

	wg.Done()
}

func ReadGenomeGraph(filename string) *SimpleGraph{
	reader := simpleio.NewSimpleReader(filename)
	defer reader.Close()
	genome := NewGraph()
	var seqIdx int32 = -1
	var i int
	var line []byte
	var done, ok bool
	var words, text [][]byte
	var weight float32
	var currNode *Node
	//var currNode *simpleGraph.Node
	//Uses this map to add edges to graph
	edges := make(map[string]*Node)


	for line, done = simpleio.ReadLine(reader); !done; line, done = simpleio.ReadLine(reader) {
		switch true {
		case bytes.HasPrefix(line, []byte(">")):
			seqIdx++
			words = bytes.Split(line[1:], []byte(":"))

			currNode = &Node{Id: uint32(seqIdx), Name: string(words[0])}
			if len(words) == 2 {
				text = bytes.Split(words[1], []byte("_"))
				currNode.Info = &Annotation{Allele: uint8(common.StringToUint32(string(text[1]))), Start: common.StringToUint32(string(text[3])), Variant: uint8(common.StringToUint32(string(text[2])))}
			}
			AddNode(genome, currNode)
			_, ok = edges[string(line[1:])]

			if !ok {
				edges[string(line[1:])] = currNode
			}
		case bytes.Contains(line, []byte("\t")):
			words = bytes.Split(line, []byte("\t"))
			if len(words) > 2 {
				for i = 1; i < len(words); i += 2 {
					weight = float32(common.StringToFloat64(string(words[i])))
					AddEdge(edges[string(words[0])], edges[string(words[i+1])], weight)
				}
			}
		case !bytes.ContainsAny(line, "\t:"):
			 genome.Nodes[seqIdx].Seq = append(genome.Nodes[seqIdx].Seq, simpleio.ByteSliceToDnaBases(line)...)

			 ///for line, done = 
		}
	}
	
	for i = 0; i < len(genome.Nodes); i++ {
		genome.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(genome.Nodes[i].Seq)
	}

	return genome
}

func percentOfFloat(part int, total int) float64 {

	return (float64(part) * float64(100)) / float64(total)
}
