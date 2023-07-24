// Package hic is written to process output from the aidenlab functions
// straw function: https://github.com/aidenlab/straw
// note: this package only reads and does not write this file type
package hic

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
)

// Straw stores the output from juicer tools straw command, start of bin 1, start of bin 2 and the contacts between them, to convert to bedpe
type Straw struct {
	Bin1Start    int
	Bin2Start    int
	contactScore int
}

// Read returns a slice of straw structs from the straw input file
func Read(filename string) []Straw {
	var line string
	var answer []Straw
	var err error
	var doneReading bool

	file := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextLine(file) {
		current := processStrawLine(line)
		answer = append(answer, current)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// processStrawLine is a helper function to read a straw file into a straw struct
func processStrawLine(line string) Straw {
	words := strings.Split(line, "\t")
	startBin1 := parse.StringToInt(words[1])
	startBin2 := parse.StringToInt(words[2])
	score := parse.StringToInt(words[3])

	current := Straw{Bin1Start: startBin1, Bin2Start: startBin2, contactScore: score}

	return current
}
