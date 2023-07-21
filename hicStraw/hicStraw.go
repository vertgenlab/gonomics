// Package hicStraw is written to process output from the aidenlab straw function: https://github.com/aidenlab/straw
// note: this package only reads and does not write this file type
package hicStraw

type Straw struct {
	Bin1Start    int
	Bin2Start    int
	contactScore int
}
