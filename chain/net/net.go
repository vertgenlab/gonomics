package net

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
)

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

func Read(filename string) []Net {
	var currTName string
	var done bool
	var block Net
	file := fileio.EasyOpen(filename)
	var answer []Net
	defer file.Close()
	for block, currTName, done = NextNet(file, currTName); !done; block, currTName, done = NextNet(file, currTName) {
		if block.Class != "" {
			answer = append(answer, block)
		}
	}
	return answer
}

func NextNet(reader *fileio.EasyReader, currTName string) (Net, string, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Net{}, currTName, true
	}
	n, currTName := NewNet(line, currTName)
	return n, currTName, false
}

func NewNet(text string, currTName string) (Net, string) {
	data := strings.Split(text, " ")
	if data[0] == "net" {
		currTName = data[1]
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
