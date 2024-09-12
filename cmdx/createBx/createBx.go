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

func createBx(length int, num int, dist int, outfile string) {
	var bx []dna.Base
	var i int
	var r float64 = 1.0 / float64(length)
	var pass bool
	bx = initializeBx(length)
	var prevBx []dna.Base = make([]dna.Base, length)
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
		modifyBase(r, bx)
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

func modifyBase(r float64, bx []dna.Base) {
	idx := choseIdx(r)
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

func choseIdx(r float64) int {
	var i int = 1
	f := rand.Float64()
	for true {
		if f <= r*float64(i) {
			return i - 1
		}
		i++
	}
	return -1
}

func choseBase() int {
	randFloat := rand.Float64()
	if randFloat < 0.33 {
		return 1
	}
	if randFloat < 0.66 {
		return 2
	}
	return 3
}

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

	createBx(*L, *N, *D, flag.Arg(0))

}
