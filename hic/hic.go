// Package hic is written to process output from the aidenlab functions
// straw function: https://github.com/aidenlab/straw
// note: this package only reads and does not write this file type
package hic

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"strings"
)

// Straw stores the output from juicer tools straw command, start of bin 1, start of bin 2 and the contacts between them, to convert to bedpe
type Straw struct {
	Bin1Start    int
	Bin2Start    int
	ContactScore int
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

	current := Straw{Bin1Start: startBin1, Bin2Start: startBin2, ContactScore: score}

	return current
}

// Equal compares two Straw structs to see if the values are identical and returns a bool
func Equal(a Straw, b Straw) bool {
	var equal bool
	if a.Bin1Start == b.Bin1Start && a.Bin2Start == b.Bin2Start && a.ContactScore == b.ContactScore {
		equal = true
	} else if a.Bin1Start == b.Bin2Start && a.Bin2Start == b.Bin1Start && a.ContactScore == b.ContactScore { //bins don't have restrictions of ordering
		equal = true
	} else {
		equal = false
	}

	return equal
}

// AllAreEqual compares two slices of Straw structs to see if the values are identical and returns a bool
func AllAreEqual(a []Straw, b []Straw) bool {
	var equal bool
	var i int

	if len(a) != len(b) {
		log.Panic("Slices of Straws are not of equal length, compare records with Equal instead of AllAreEqual.")
	}

	for i = range a {
		equal = Equal(a[i], b[i])
		if !equal {
			return equal
		}
	}

	return equal
}
