package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"math/rand"
	"sort"
)

// createBx is a function that selects dna barcodes. The number of barcodes, the length (bp) and how many bases must differ
// between barcodes are all user-input variables. It works by initializing a barcode, it then randomly mutates the sequence
// and tests if it has created an acceptable barcode every few iterations. Homopolymers >= 3 are not accepted
func createBx(length int, num int, dist int, outfile string) {
	var bx, prevBx []dna.Base = initializeBx(length), make([]dna.Base, length)
	var i int
	var proportion float64 = 1.0 / float64(length)
	var pass bool
	bx = initializeBx(length)
	mp := map[string]int{
		dna.BasesToString(bx): 0,
	}

	for true {
		i++
		if i%dist == 0 {
			pass = testBx(mp, bx, dist)
			if pass {
				mp[dna.BasesToString(bx)] = 0
			}
		}
		copy(prevBx, bx)
		modifyBase(proportion, bx)
		if dna.TestForHomopolymer(bx, 3) {
			copy(bx, prevBx)
			i--
		}
		if len(mp) >= num {
			break
		}
	}
	writeRes(mp, outfile)
}

// this function writes the output in fasta format
func writeRes(mp map[string]int, outfile string) {
	var writeSlice []string
	o := fileio.EasyCreate(outfile)
	for i := range mp {
		writeSlice = append(writeSlice, i)
	}
	sort.Strings(writeSlice)
	for i := range writeSlice {
		fileio.WriteToFileHandle(o, fmt.Sprintf(">bx%d\n%s", i, writeSlice[i]))
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

// testBx tests to see if the current barcode is an acceptable barcode to add to the output. It takes the map of currently
// accepted barcodes, the test barcodes, and the user variable of how many bases must differ between barcodes
func testBx(mp map[string]int, bx []dna.Base, dist int) bool {
	_, found := mp[dna.BasesToString(bx)]
	if found {
		return false
	}
	for i := range mp {
		if dna.Dist(dna.StringToBases(i), bx) < dist {
			return false
		}
	}
	return true
}

// modifyBase is the function which alters the barcode from one iteration of the loop to the next. It has helper functions
// that select the index and the new identity of a single base. It edits the bx slice directly
func modifyBase(proportion float64, bx []dna.Base) {
	idx := choseIdx(proportion)
	mod := choseBase()
	switch int(bx[idx]) + mod {
	case 1, 5:
		bx[idx] = 1
	case 2, 6:
		bx[idx] = 2
	case 3:
		bx[idx] = 3
	case 4:
		bx[idx] = 0
	}
}

// choseIdx uses a random float to pick the base at what index of the barcode that will be modified. It takes a "proportion"
// which is calculated with 1.0 / length(barcode). The index (i) is selected with: i * proportion < randomFloat < i+1 * proportion.
func choseIdx(proportion float64) int {
	var idx int = 1
	randfloat := rand.Float64()
	for true {
		if randfloat <= proportion*float64(idx) {
			return idx - 1
		}
		idx++
	}
	return -1
}

// choseBase uses a random float to pick if the base to be modified will be iterated by 1, 2, or 3
func choseBase() int {
	randfloat := rand.Float64()
	if randfloat <= 0.33 {
		return 1
	}
	if randfloat <= 0.66 {
		return 2
	}
	return 3
}

// initializeBx creates the initial barcode by repeating ACGT... until the length is satisfied
func initializeBx(length int) []dna.Base {
	var bx []dna.Base
	var b dna.Base = 0
	for i := 0; i < length; i++ {
		bx = append(bx, b)
		if b == dna.StringToBase("T") {
			b = 0
		} else {
			b++
		}
	}
	return bx
}

func main() {
	var L *int = flag.Int("L", 5, "Length of barcode to create")
	var N *int = flag.Int("N", 10, "Number of barcodes to create")
	var D *int = flag.Int("D", 2, "Number of bases which must be different between acceptable barcodes")
	flag.Parse()

	if len(flag.Args()) != 1 {
		log.Fatalf("Expected 1 argument, got %d", len(flag.Args()))
	}
	outfile := flag.Arg(0)
	createBx(*L, *N, *D, outfile)

}
