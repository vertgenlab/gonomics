package chromInfo

import (
	"bufio"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type ChromInfo struct {
	Name  string
	Size  int64
	Order int64
}

func SliceToMap(chroms []*ChromInfo) map[string]*ChromInfo {
	answer := make(map[string]*ChromInfo)
	for i := 0; i < len(chroms); i++ {
		curr := ChromInfo{Name: chroms[i].Name, Size: chroms[i].Size, Order: chroms[i].Order}
		answer[chroms[i].Name] = &curr
	}
	return answer
}

func ReadToSlice(filename string) []*ChromInfo {
	var line string
	var answer []*ChromInfo
	var count int64
	var doneReading bool

	file := fileio.MustOpen(filename)
	defer file.Close()
	reader := bufio.NewReader(file)

	count = 0
	for line, doneReading = fileio.NextRealLine(reader); !doneReading; line, doneReading = fileio.NextRealLine(reader) {
		words := strings.Fields(line)
		if len(words) != 2 {
			log.Fatalf("Error: expecting 2 columns, but got %d on line:%s\n", len(words), line)
		}
		size := common.StringToInt64(words[1])
		curr := ChromInfo{Name: words[0], Size: size, Order: count}
		answer = append(answer, &curr)
		count++
	}
	return answer
}

func ReadToMap(filename string) map[string]*ChromInfo {
	var line string
	answer := make(map[string]*ChromInfo)
	var count int64
	var doneReading bool

	file := fileio.MustOpen(filename)
	defer file.Close()
	reader := bufio.NewReader(file)

	count = 0
	for line, doneReading = fileio.NextRealLine(reader); !doneReading; line, doneReading = fileio.NextRealLine(reader) {
		words := strings.Fields(line)
		if len(words) != 2 {
			log.Fatalf("Error: expecting 2 columns, but got %d on line:%s\n", len(words), line)
		}
		size := common.StringToInt64(words[1])
		curr := ChromInfo{Name: words[0], Size: size, Order: count}
		answer[words[0]] = &curr
		count++
	}
	return answer
}
