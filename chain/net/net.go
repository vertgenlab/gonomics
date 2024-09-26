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
}

// Read parses a net file and returns a slice of Net and a map with chromosome names and sizes in the net file
func Read(filename string) ([]Net, map[string]chromInfo.ChromInfo) {
	var currTName string
	var done bool
	var block Net
	mp := make(map[string]chromInfo.ChromInfo)
	file := fileio.EasyOpen(filename)
	var answer []Net
	defer file.Close()
	for block, currTName, done = NextNet(file, currTName, mp); !done; block, currTName, done = NextNet(file, currTName, mp) {
		if block.Class != "" {
			answer = append(answer, block)
		}
	}
	return answer, mp
}

// NextNet is a helper function for the Read function. It reads lines from an EasyReader and returns a Net struct, the
// current chromosome, and a bool which will return true once it has reached the end of the file
func NextNet(reader *fileio.EasyReader, currTName string, mp map[string]chromInfo.ChromInfo) (Net, string, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Net{}, currTName, true
	}
	n, currTName := NewNet(line, currTName, mp)
	return n, currTName, false
}

// NewNet is a helper function for the Read function. It takes in a single line of the file representing a Net entry
// and returns a filled in Net struct and the current chromosome
func NewNet(text string, currTName string, mp map[string]chromInfo.ChromInfo) (Net, string) {
	data := strings.Split(text, " ")
	if data[0] == "net" {
		currTName = data[1]
		mp[data[1]] = chromInfo.ChromInfo{Name: data[1], Size: parse.StringToInt(data[2])}
		return Net{}, currTName
	}
	lvl, newData := countLevel(data)
	curr := Net{
		TName:       currTName,
		Level:       lvl,
		Class:       newData[0],
		TStart:      parse.StringToInt(newData[1]),
		TSize:       parse.StringToInt(newData[2]),
		QName:       newData[3],
		Orientation: parse.StringToStrand(newData[4]),
		QStart:      parse.StringToInt(newData[5]),
		QSize:       parse.StringToInt(newData[6]),
		ExtraFields: strings.Join(newData[7:], " "),
	}

	if curr.Class == "gap" {
		curr.Level--
	}
	if curr.Level > 1 {
		curr.Level--
	}
	return curr, currTName
}

// countLevel is a helper struct for the NewNet function. It parses how many spaces proceed the data on the line
// of the file to determine the Level of the net. It also returns the data without the proceeding spaces for easy parsing
func countLevel(data []string) (int, []string) {
	var c int
	for i := range data {
		if data[i] == "" {
			c++
		} else {
			return c, data[i:]
		}
	}
	return -1, []string{}
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
		fileio.WriteToFileHandle(file, ToString(nets[i]))
		prevChrom = currChrom
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// ToString takes a Net struct and formats it for writing
func ToString(n Net) string {
	var rec string = fmt.Sprintf("%s %d %d %s %c %d %d %s", n.Class, n.TStart, n.TSize, n.QName, parse.StrandToRune(n.Orientation), n.QStart, n.QSize, n.ExtraFields)
	for i := 0; i < n.Level; i++ {
		rec = fmt.Sprintf(" %s", rec)
	}
	if n.Class == "gap" {
		rec = fmt.Sprintf(" %s", rec)
	}
	if n.Level > 1 {
		rec = fmt.Sprintf(" %s", rec)
	}
	return rec
}

func GoReadToChan(filename string) (<-chan Net, map[string]chromInfo.ChromInfo) {
	var wg sync.WaitGroup
	var currTName string
	chromSizes := make(map[string]chromInfo.ChromInfo)
	file := fileio.EasyOpen(filename)
	data := make(chan Net, 1000)
	wg.Add(1)
	go ReadToChan(file, data, currTName, chromSizes, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, chromSizes
}

func ReadToChan(file *fileio.EasyReader, data chan<- Net, currTName string, chromSizes map[string]chromInfo.ChromInfo, wg *sync.WaitGroup) {
	for curr, currTName, done := NextNet(file, currTName, chromSizes); !done; curr, currTName, done = NextNet(file, currTName, chromSizes) {
		if curr.Class != "" {
			data <- curr
		}
	}
	file.Close()
	wg.Done()
}
