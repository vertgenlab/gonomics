package qDna

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"math"
	"os"
	"strings"
)

type ChrDict struct {
	Chr   string
	Coord int64
	//Size int64
}

func calcMask(n int) int64 {
	var answer float64 = 0
	for i := 0; i < n; i++ {
		answer += math.Exp2(float64(i))
	}
	return int64(answer)
}

func putTogether(seq []dna.Base, start int, end int) int64 {
	var answer int64 = int64(seq[start])
	for i := start + 1; i < end; i++ {
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

func IndexRefSlidingWindow(ref []*QFrag, seed int) map[int64][]int64 {
	answer := make(map[int64][]int64)
	var i, nodeIdx int64

	var sequence []dna.Base
	//var curr ChrDict
	for i = 0; i < int64(len(ref)); i++ {
		sequence = mostLikelySeq(ref[i].Seq)
		nodeIdx = i << 32
		for j := 0; j < len(sequence)-seed+1; j++ {
			//curr = ChrDict{Chr: ref[i].From[0].Chr, Coord: int64(j)}
			answer[putTogether(sequence, j, j+seed)] = append(answer[putTogether(sequence, j, j+seed)], nodeIdx|int64(j))
		}
	}

	return answer
}

func ReturnIdPos(code int64) (int64, int64) {
	//2^33-1 = 8589934591 or 0000000000000000000000000000000011111111111111111111111111111111
	var leftMask int64 = 8589934591
	//var rightMask int64 =1111111111111111111111111111111100000000000000000000000000000000
	var rightMask int64 = leftMask << 32
	var nodeId, pos int64
	nodeId = code | rightMask
	nodeId = nodeId >> 32
	pos = code | leftMask
	return nodeId, pos
}

func Read(filename string) map[int64][]int64 {
	answer := make(map[int64][]int64)
	file, _ := os.Open(filename)
	defer file.Close()
	reader := bufio.NewReader(file)
	var err error
	var line string

	var words []byte
	for ; err != io.EOF; words, _, err = reader.ReadLine() {
		line = string(words[:])
		data := strings.Split(line, " ")
		for i := 1; i < len(data)-1; i++ {
			answer[common.StringToInt64(data[0])] = append(answer[common.StringToInt64(data[0])], common.StringToInt64(data[i+1]))
		}
		//fmt.Println(data)
	}
	return answer
}

func WriteDictToFileHandle(file *os.File, input map[int64][]int64) error {
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

func Write(filename string, data map[int64][]int64) {
	file := fileio.MustCreate(filename)
	defer file.Close()
	WriteDictToFileHandle(file, data)
}
