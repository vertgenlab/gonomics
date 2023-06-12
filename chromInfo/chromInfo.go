// Package chromInfo provides function to read and manipulate chromInfo files
package chromInfo

import (
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// ChromInfo stores the chromosome name and size of each chromosome as well as the order they appear in the file.
type ChromInfo struct {
	Name  string
	Size  int
	Order int
}

// SliceToMap converts a slice of chromInfo structs into a map[Name]ChromInfo.
func SliceToMap(chroms []ChromInfo) map[string]ChromInfo {
	answer := make(map[string]ChromInfo)
	for i := 0; i < len(chroms); i++ {
		curr := ChromInfo{Name: chroms[i].Name, Size: chroms[i].Size, Order: chroms[i].Order}
		answer[chroms[i].Name] = curr
	}
	return answer
}

// ReadToSlice reads a chromInfo file into a slice where the index corresponds to the order of appearance in the file.
func ReadToSlice(filename string) []ChromInfo {
	var line string
	var answer []ChromInfo
	var count int
	var doneReading bool

	file := fileio.EasyOpen(filename)

	count = 0
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		answer = append(answer, readLine(line, count))
		count++
	}

	err := file.Close()
	exception.PanicOnErr(err)
	return answer
}

// SliceToMap reads a chromInfo file into a map[Name]ChromInfo.
func ReadToMap(filename string) map[string]ChromInfo {
	data := ReadToSlice(filename)
	return SliceToMap(data)
}

// readLine parses a single line of a chromInfo file and returns a ChromInfo struct.
func readLine(line string, count int) ChromInfo {
	var answer ChromInfo

	words := strings.Fields(line)
	if len(words) != 2 {
		log.Panicf("Error: expecting 2 columns, but got %d on line:%s\n", len(words), line)
	}

	size, err := strconv.Atoi(words[1])
	exception.PanicOnErr(err)

	answer = ChromInfo{Name: words[0], Size: size, Order: count}
	return answer
}
