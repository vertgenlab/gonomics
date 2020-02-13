package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"math/rand"
)

func GenomeDiversity(genome []*Node, readLength int, numReads int, numChanges int, snpInsDel bool) []*fastq.PairedEnd {
	var answer []*fastq.PairedEnd = make([]*fastq.PairedEnd, numReads)
	var start int
	var chromIdx int
	var strand bool
	var editGenome int
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateDiverseFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			if snpInsDel {
				editGenome = randIntInRange(0, 8)
				//snp
				if editGenome == 0 || editGenome == 1 {
					mutate(curr.Seq, numChanges)
				} else if editGenome == 2 {
					curr = deletions(curr, randIntInRange(0, 3))
				} else if editGenome == 3 || editGenome == 4 {
					curr = insertion(curr, randIntInRange(0, 3))
				} else {
					//no mutations
				}
			}
			pairEnd := fastq.PairedEnd{}
			pairEnd = singleToPair(genome, curr)

			answer[i] = &pairEnd
			i++
		}
	}
	return answer
}

func singleToPair(genome []*Node, inFq fastq.Fastq) fastq.PairedEnd {
	fqPair := fastq.PairedEnd{Fwd: nil, Rev: nil}
	fqPair.Fwd = fastq.Copy(&inFq)
	fqPair.Rev = fastq.Copy(fqPair.Fwd)
	fastq.ReverseComplement(fqPair.Rev)
	fqPair.Fwd.Name = fmt.Sprintf("%s_R1", fqPair.Fwd.Name)
	fqPair.Rev.Name = fmt.Sprintf("%s_R2", fqPair.Rev.Name)
	return fqPair
}

/*
func GenomeDiversity(genome []*Node, readLength int, numOfReads int, edits int) []*fastq.Fastq {
	//weights for ranomization
	var snp int = 8*edits
	var ins int = 2*edits
	var del int = 4*edits
	var mutations int
	var gDna []*fastq.Fastq = make([]*fastq.Fastq, numOfReads)

	for monsters := 0; monsters < numOfReads; {
		mutations = randIntInRange(0, snp+ins+del+8)
		//ins, del, snp
		if mutations < 0 && mutations < ins {
			gDna[monsters] = insertions(RandomReadOneMutation(genome, readLength, 0), ins)
		} else if mutations < ins+del && mutations > del {
			gDna[monsters] = deletions(RandomReadOneMutation(genome, readLength+del, 0), del)
		} else if mutations < snp+del+del && mutations > del+ins {
			gDna[monsters] = RandomReadOneMutation(genome, readLength, 2)
		} else {
			gDna[monsters] = RandomReadOneMutation(genome, readLength, 0)
		}
	}
	return gDna
}*/

func RandomReads(genome []*Node, readLength int, numReads int, numChanges int, snpInsDel bool) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var strand bool
	var editGenome int
	//size of insertion or deletion
	var size int = 5
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength-size)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			if snpInsDel {
				editGenome = randIntInRange(0, 30)
				//del
				if editGenome < 4 {
					location := int64(randIntInRange(size, len(curr.Seq)-size))
					dna.Delete(curr.Seq, location, location+int64(size))
					curr.Name = fmt.Sprintf("%s_del:%dto%d", curr.Name, location, location+int64(size))
					curr.Seq = append(curr.Seq, genome[chromIdx].Seq[start+readLength:size+start+readLength]...)
					//log.Printf("Size of read with deletions is=%d\n", len(curr.Seq))

				} else if editGenome < 8 && editGenome > 4 {
					curr = insertion(curr, 2)
				} else if editGenome < 8 && editGenome < 24 {
					mutate(curr.Seq, numChanges)
					mutate(curr.Seq, numChanges)
				} else {
					//no mutations
				}

			}
			if len(curr.Seq) != readLength {
				log.Fatalf("incoreect read len %d", len(curr.Seq))
			}
			answer[i] = &curr
			i++

		}
	}

	return answer
}

func deletions(del fastq.Fastq, size int) fastq.Fastq {
	var location int = randIntInRange(size, len(del.Seq)-size-1)
	undeleted := make([]dna.Base, len(del.Seq)-size)
	copy(undeleted[0:location], del.Seq[location+size:len(del.Seq)])
	del.Name = fmt.Sprintf("%s_del:%dto%d", del.Name, location, location+size)
	del.Seq = undeleted
	return del
}

func insertion(ins fastq.Fastq, size int) fastq.Fastq {
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	extra := make([]dna.Base, size)
	for base := 0; base < size; base++ {
		extra = append(extra, possibleBases[randIntInRange(0, len(possibleBases))])
	}
	var location int = randIntInRange(size, len(ins.Seq)-size-1)
	copy(ins.Seq[location:location+len(extra)], extra[0:len(extra)])
	ins.Name = fmt.Sprintf("%s_ins:%dto%d", ins.Name, location, location+size)
	return ins
}

func RandomFastqGen(genome []*fasta.Fasta, readLength int, numReads int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var readName string
	var strand bool
	var qual []rune
	var seq []dna.Base
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			readName = fmt.Sprintf("%s_%d_%d_%c", genome[chromIdx].Name, start, start+readLength, common.StrandToRune(strand))
			qual = generateDiverseFakeQual(readLength)
			seq = genome[chromIdx].Seq[start : start+readLength]
			if !strand {
				dna.ReverseComplement(seq)
			}
			answer[i] = &fastq.Fastq{Name: readName, Seq: seq, Qual: qual}
			i++
		}
	}
	return answer
}

func randIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func mutate(sequence []dna.Base, numChanges int) {
	possibleBases := []dna.Base{0, 1, 2, 3}
	for i := 0; i < numChanges; i++ {
		sequence[randIntInRange(0, len(sequence))] = possibleBases[randIntInRange(0, len(possibleBases))]
	}
}

func mutatePos(seq []dna.Base, pos int) {
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	newBase := dna.A
	for newBase = possibleBases[randIntInRange(0, len(possibleBases))]; newBase != seq[pos]; newBase = possibleBases[randIntInRange(0, len(possibleBases))] {
	}
	seq[pos] = newBase
}

func RandomReadOneMutation(genome []*Node, readLength int, mutantPos int) *fastq.Fastq {
	var answer *fastq.Fastq = nil
	var start int
	var chromIdx int
	var strand bool
	var readOk bool = false

	for !readOk {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			strand = randIntInRange(0, 2) == 0
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			mutatePos(curr.Seq, mutantPos)
			answer = &curr
			readOk = true
		}
	}
	return answer
}

func generateDiverseFakeQual(length int) []rune {
	var answer []rune = make([]rune, length)
	var asci = []rune{'!', '#', '$', '%', '&', '(', ')', '*', '+', '`', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}

	for i := 0; i < length; i++ {
		answer[i] = asci[randIntInRange(0, len(asci))]
	}
	return answer
}

func generateFakeQual(length int) []rune {
	var answer []rune = make([]rune, length)
	for i := 0; i < length; i++ {
		answer[i] = 'J'
	}
	return answer
}

func generateQual(bases []dna.Base) []rune {
	var ans []rune
	for i := 0; i < len(bases); i++ {
		ans = append(ans, 'J')
	}
	return ans
}
