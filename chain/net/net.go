package net

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
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

func NextNet(reader *fileio.EasyReader, currTName string, mp map[string]chromInfo.ChromInfo) (Net, string, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Net{}, currTName, true
	}
	n, currTName := NewNet(line, currTName, mp)
	return n, currTName, false
}

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
