package net

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
	"sync"
)

// Net represents an alignment block built off of chains, it represents the UCSC net file type
type Net struct {
	TName       string
	Level       int
	Class       string
	TStart      int
	TSize       int
	QName       string
	Orientation bool
	QStart      int
	QSize       int
	ExtraFields string
	spaces      int //used for writing
}

// Read parses a net file and returns a slice of Net and a map with chromosome names and sizes in the net file
func Read(filename string) ([]Net, map[string]chromInfo.ChromInfo) {
	var currTName string
	var done bool
	var block Net
	var currLvl [][]int = [][]int{{0, 0}, {0}} //{level, spaces}, {lvlKey}
	mp := make(map[string]chromInfo.ChromInfo)
	file := fileio.EasyOpen(filename)
	var answer []Net
	defer file.Close()
	for block, currTName, currLvl, done = NextNet(file, currTName, currLvl, mp); !done; block, currTName, currLvl, done = NextNet(file, currTName, currLvl, mp) {
		if block.Class != "" {
			answer = append(answer, block)
		}
	}
	return answer, mp
}

// NextNet is a helper function for the Read function. It reads lines from an EasyReader and returns a Net struct, the
// current chromosome, and a bool which will return true once it has reached the end of the file
func NextNet(reader *fileio.EasyReader, currTName string, currLvl [][]int, mp map[string]chromInfo.ChromInfo) (Net, string, [][]int, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Net{}, currTName, currLvl, true
	}
	n, currTName, currLvl := NewNet(line, currTName, currLvl, mp)
	return n, currTName, currLvl, false
}

// NewNet is a helper function for the Read function. It takes in a single line of the file representing a Net entry
// and returns a filled in Net struct and the current chromosome
func NewNet(text string, currTName string, currLvl [][]int, mp map[string]chromInfo.ChromInfo) (Net, string, [][]int) {
	data := strings.Split(text, " ")
	if data[0] == "net" {
		currTName = data[1]
		mp[data[1]] = chromInfo.ChromInfo{Name: data[1], Size: parse.StringToInt(data[2])}
		return Net{}, currTName, [][]int{{0, 0}, {0}}
	}
	currLvl, newData := determineLevel(data, currLvl)
	curr := Net{
		TName:       currTName,
		Level:       currLvl[0][0],
		Class:       newData[0],
		TStart:      parse.StringToInt(newData[1]),
		TSize:       parse.StringToInt(newData[2]),
		QName:       newData[3],
		Orientation: parse.StringToStrand(newData[4]),
		QStart:      parse.StringToInt(newData[5]),
		QSize:       parse.StringToInt(newData[6]),
		ExtraFields: strings.Join(newData[7:], " "),
		spaces:      currLvl[0][1],
	}

	return curr, currTName, currLvl
}

// determineLevel is a helper struct for the NewNet function. It parses how many spaces proceed the data on the line
// of the file to determine the Level of the net. It also returns the data without the proceeding spaces for easy parsing
func determineLevel(data []string, currLvl [][]int) ([][]int, []string) {
	var spaces int
	for i := range data {
		if data[i] == "" {
			spaces++
		} else {
			break
		}
	}
	data = data[spaces:]

	if spaces == currLvl[0][1] {
		return currLvl, data
	} else if spaces > currLvl[0][1] {
		if data[0] == "fill" {
			currLvl[0][0]++
		}
		currLvl[1] = append(currLvl[1], currLvl[0][0])
		currLvl[0][1] = spaces
		return currLvl, data
	} else if spaces < currLvl[0][1] {
		currLvl[0][1] = spaces
		if data[0] == "fill" {
			currLvl[0][0] = currLvl[1][spaces]
		} else if data[0] == "gap" {
			currLvl[0][0] = currLvl[1][spaces-1]
		}
		currLvl[1] = currLvl[1][:spaces]
		return currLvl, data
	}

	return [][]int{}, []string{}
}

// Write takes in a file name, a slice of Net and a map of chromInfo and writes out a properly formatted net file.
func Write(outfile string, nets []Net, chromSizes map[string]chromInfo.ChromInfo) {
	var currChrom, prevChrom string
	file := fileio.EasyCreate(outfile)
	for i := range nets {
		currChrom = nets[i].TName
		if currChrom != prevChrom {
			fileio.WriteToFileHandle(file, fmt.Sprintf("net %s %d", currChrom, chromSizes[currChrom].Size))
		}
		fileio.WriteToFileHandle(file, ToString(nets[i], true))
		prevChrom = currChrom
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// ToString takes a Net struct and formats it for writing. If spaces is set to true, the proceeding spaces will be added according to the spaces field in the Net struct
func ToString(n Net, spaces bool) string {
	var sb strings.Builder
	if !spaces {
		return fmt.Sprintf("%s %d %d %s %c %d %d %s", n.Class, n.TStart, n.TSize, n.QName, parse.StrandToRune(n.Orientation), n.QStart, n.QSize, n.ExtraFields)
	}
	str := fmt.Sprintf("%s %d %d %s %c %d %d %s", n.Class, n.TStart, n.TSize, n.QName, parse.StrandToRune(n.Orientation), n.QStart, n.QSize, n.ExtraFields)
	for i := 0; i < n.spaces; i++ {
		sb.WriteString(" ")
	}
	sb.WriteString(str)
	return sb.String()
}

func GoReadToChan(filename string) (<-chan Net, map[string]chromInfo.ChromInfo) {
	var wg sync.WaitGroup
	var currTName string
	var currLvl [][]int = [][]int{{0, 0}, {0}}
	chromSizes := make(map[string]chromInfo.ChromInfo)
	file := fileio.EasyOpen(filename)
	data := make(chan Net, 1000)
	wg.Add(1)
	go ReadToChan(file, data, currTName, currLvl, chromSizes, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, chromSizes
}

func ReadToChan(file *fileio.EasyReader, data chan<- Net, currTName string, currLvl [][]int, chromSizes map[string]chromInfo.ChromInfo, wg *sync.WaitGroup) {
	for curr, currTName, currLvl, done := NextNet(file, currTName, currLvl, chromSizes); !done; curr, currTName, currLvl, done = NextNet(file, currTName, currLvl, chromSizes) {
		if curr.Class != "" {
			data <- curr
		}
	}
	file.Close()
	wg.Done()
}
