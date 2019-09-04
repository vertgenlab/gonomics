package qDna

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
	"bufio"
	"math"
	"os"
)

type ChrDict struct {
	Chr   string
	Coord int64
	//Size int64
}

//Creates a hash table using encoded bytes seed must be less than 32
func indexChr(chr *QFrag, seed int) map[int64][]int64 {
	answer := make(map[int64][]int64)

	var sequence []dna.Base = mostLikelySeq(chr.Seq)
	dna.AllToUpper(sequence)
	//seedLength, tail := int64(len(sequence))/denominator, numerator%denominator
	leftover := len(chr.Seq) % seed
	for i := 0; i < len(sequence)-leftover; i += seed {
		answer[putTogether(sequence[i:i+seed])] = append(answer[putTogether(sequence[i:i+seed])], int64(i))
	}
	if leftover != 0 {
		answer[putTogether(sequence[len(chr.Seq)-leftover:len(chr.Seq)])] = append(answer[putTogether(sequence[len(chr.Seq)-leftover:len(chr.Seq)])], int64(len(chr.Seq)-leftover))
	}
	return answer
}

func indexRef(ref []*QFrag, seed int) map[int64][]ChrDict {
	answer := make(map[int64][]ChrDict)

	var sequence []dna.Base
	var leftover int
	for i := 0; i < len(ref); i++ {

		sequence = mostLikelySeq(ref[i].Seq)
		leftover = len(ref[i].Seq) % seed
		for j := 0; j < len(sequence)-leftover; j += seed {
			answer[putTogether(sequence[j:j+seed])] = append(answer[putTogether(sequence[j:j+seed])], ChrDict{Chr: ref[i].From[0].Chr, Coord: int64(j)})
		}
		if leftover != 0 {
			//logic to deal with the tail
			answer[putTogether(sequence[len(sequence)-leftover:len(sequence)])] = append(answer[putTogether(sequence[len(sequence)-leftover:len(sequence)])], ChrDict{Chr: ref[i].From[0].Chr, Coord: int64(len(sequence) - leftover)})
		}
	}
	return answer
}

func calcMask(n int) int64 {
	var answer float64 = 0
	for i := 0; i < n; i++ {
		answer += math.Exp2(float64(i))
	}
	return int64(answer)
}

func putTogether(seq []dna.Base) int64 {
	var answer int64 = int64(seq[0])
	for i := 1; i < len(seq); i++ {
		//fmt.Printf("answer is (before): %b\n", answer)
		answer = answer << 2
		//fmt.Printf("answer is (halfway): %b\n", answer)
		answer = answer | int64(seq[i])
		//fmt.Printf("answer is (loopEnd): %b\n", answer)
	}
	return answer
}

func ChromMap(ref []*QFrag) map[string][]*QBase {
	m := make(map[string][]*QBase)
	//var answer []*fasta.Fasta
	var curr *QFrag
	for i := 0; i < len(ref); i++ {
		curr = ref[i]
		_, ok := m[curr.From[0].Chr]
		if !ok {
			m[curr.From[0].Chr] = curr.Seq
		}
	}
	return m
}

func IndexRefSlidingWindow(ref []*QFrag, seed int) map[int64][]ChrDict {
	answer := make(map[int64][]ChrDict)
	var sequence []dna.Base
	//var curr ChrDict
	for i := 0; i < len(ref); i++ {
		sequence = mostLikelySeq(ref[i].Seq)
		for j := 0; j < len(sequence)-seed; j++ {
			//curr = ChrDict{Chr: ref[i].From[0].Chr, Coord: int64(j)}
			answer[putTogether(sequence[j:j+seed])] = append(answer[putTogether(sequence[j:j+seed])], ChrDict{Chr: ref[i].From[0].Chr, Coord: int64(j)})
		}
	}
	return answer
}

func Read(filename string) map[int64][]ChrDict {
	answer := make(map[int64][]ChrDict)
	file, _ := os.Open(filename)
	defer file.Close()
	reader := bufio.NewReader(file)

	var err error
	var line string
	//var currNode *Node
	var words []byte
	for ; err != io.EOF; words, _, err = reader.ReadLine() {
		line = string(words[:])
		data := strings.Split(line, " ")
		for i := 1; i < len(data)-1; i++ {
			answer[common.StringToInt64(data[0])] = append(answer[common.StringToInt64(data[0])], ChrDict{Chr: data[i], Coord: common.StringToInt64(data[i+1])})
		}
		fmt.Println(data)

	}
	return answer
}

func WriteDictToFileHandle(file *os.File, input map[int64][]ChrDict) error {
	var err error

	for i := range input {
		_, err = fmt.Fprintf(file, "%v", i)
		for j := range input[i] {
			_, err = fmt.Fprintf(file, "\t%v\t%s", input[i][j].Coord, input[i][j].Chr)
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
	}
	return err
}

func Write(filename string, data map[int64][]ChrDict) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	WriteDictToFileHandle(file, data)
}
