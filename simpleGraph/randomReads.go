package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"math/rand"
)

func GenomeDiversity(genome []*Node, readLength int, numReads int, numChanges int, snpInsDel bool) []*fastq.PairedEnd {
	var answer []*fastq.PairedEnd = make([]*fastq.PairedEnd, numReads)
	var start int
	var chromIdx int
	var strand bool
	var editGenome int
	for i := 0; i < numReads; {
		var size int = randIntInRange(8, 24)
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d:%d:%d:%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateDiverseFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			if snpInsDel {
				location := int64(randIntInRange(2, len(curr.Seq)-24))
				editGenome = randIntInRange(0, 30)
				//del
				if editGenome < 4 {

					dna.Delete(curr.Seq, location, location+int64(size))
					curr.Name = fmt.Sprintf("%s:del:%dto%d", curr.Name, location, location+int64(size))
					curr.Seq = append(curr.Seq, genome[chromIdx].Seq[int64(start)+location+int64(size):readLength+size]...)
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
	var size int = 0
	for i := 0; i < numReads; {
		size = randIntInRange(2, 8)
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(30, len(genome[chromIdx].Seq)-readLength-size)
		strand = randIntInRange(0, 2) == 0
		location := int64(randIntInRange(8, readLength-size))
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d:%d:%d:%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)

			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			if snpInsDel {

				editGenome = randIntInRange(0, 30)
				//del
				if editGenome <= 4 {
					//location := int64(randIntInRange(2+size, readLength-size)-8)
					deletion := fastq.Fastq{}
					deletion.Name = fmt.Sprintf("%d:%d:%d:%c:deletion:%dto%d", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand), location, location+int64(size))
					deletion.Seq = make([]dna.Base, readLength)
					copy(deletion.Seq, genome[chromIdx].Seq[start:start+readLength+size])
					//log.Printf("Sequence before deletion at location=%d\n%s\n", location, dna.BasesToString(deletion.Seq)[:10])
					dna.Delete(deletion.Seq, location, location+int64(size))
					//log.Printf("%s\nSequence after deletion len=%d\n", dna.BasesToString(deletion.Seq)[:10], len(deletion.Seq))
					//append(curr.Seq[:location], curr.Seq[location+int64(size):]...)

					curr.Name = deletion.Name
					curr.Seq = make([]dna.Base, readLength)
					copy(curr.Seq, deletion.Seq)
					//log.Printf("Sequence sending off to aligner!=%slen=%d\n", dna.BasesToString(curr.Seq), len(curr.Seq))
					///curr = deletion

					//curr.Name = fmt.Sprintf("%s_del:%dto%d", curr.Name, location, location+int64(size))

					//log.Printf("Size of read with deletions is=%d\n", len(curr.Seq))
				} else if editGenome <= 8 && editGenome > 4 {

					possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
					extra := make([]dna.Base, size)
					for base := 0; base < size; base++ {
						extra = append(extra, possibleBases[randIntInRange(0, len(possibleBases))])
					}

					insertion := fastq.Fastq{}
					insertion.Name = fmt.Sprintf("%d:%d:%d:%c:insertion:%dto%d", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand), location, location+int64(size))
					insertion.Seq = make([]dna.Base, readLength)
					copy(insertion.Seq, genome[chromIdx].Seq[start:start+readLength-size])
					//log.Printf("Sequence before insertion at location=%d\n%s\n", location, dna.BasesToString(insertion.Seq)[:10])
					dna.Insert(insertion.Seq, location, extra)
					//log.Printf("%s\nSequence after insertion len=%d\n", dna.BasesToString(insertion.Seq)[:10], len(insertion.Seq))
					//append(curr.Seq[:location], curr.Seq[location+int64(size):]...)
					curr.Name = insertion.Name
					curr.Seq = make([]dna.Base, readLength)
					copy(curr.Seq, insertion.Seq)
				} else if editGenome > 8 && editGenome <= 24 {
					snp := fastq.Fastq{}
					snp.Seq = make([]dna.Base, readLength)
					copy(snp.Seq, genome[chromIdx].Seq[start:start+readLength])
					before := snp.Seq[location]
					newMutatePos(snp.Seq, int(location))
					name := fmt.Sprintf("%d_%d_%d_SNP:%d:%c_to_%c", genome[chromIdx].Id, start, start+readLength, location, dna.BaseToRune(before), dna.BaseToRune(snp.Seq[location]))
					location = int64(randIntInRange(8+size, len(snp.Seq)-size))
					newMutatePos(snp.Seq, int(location))
					name += fmt.Sprintf("%s_&_SNP:%d:%c_to_%c_%c", snp.Name, location, dna.BaseToRune(before), dna.BaseToRune(snp.Seq[location]), common.StrandToRune(strand))
					curr.Name = name
					curr.Seq = make([]dna.Base, readLength)
					copy(curr.Seq, snp.Seq)
				} else {

				}
			}
			if len(curr.Seq) != readLength {
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

func newMutatePos(seq []dna.Base, pos int) {
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	newBase := possibleBases[randIntInRange(0, len(possibleBases))]
	if newBase != seq[pos] {
		seq[pos] = newBase
	} else {
		newMutatePos(seq, pos)
	}
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
	//var asci = []rune{'!', '#', '$', '%', '&', '(', ')', '*', '+', '`', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}
	var asci = []rune{'F', ',', 'F', ':'}
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
