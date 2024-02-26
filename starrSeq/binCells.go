package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"math/rand"
	"sort"
)

// DetermineBin takes in an int and the probablity of any given bin and returns an int corresponding to the bin the cell belongs to.
func DetermineBin(numBins int, prob float64) int {
	var bin int

	randNum := rand.Float64()      // draw a random number
	for j := 0; j < numBins; j++ { //iterate over the number of bins
		if j == numBins { //if we've gotten to the last bin and the cell still hasn't been partitioned, we don't have to check anything else we can just put it there.
			return j
		}
		if randNum >= prob*float64(j) && randNum < prob*float64(j+1) { // check to see if our random number belongs to that bin
			bin = j
		}
	}
	return bin
}

// DistributeCells takes in a slice of string that is created in parseBam that has a list of cellBarcodes and constructs. It will partition those cells and constructs into separate slices. The number of bins is a user-input variable
func DistributeCells(s ScStarrSeqSettings, cellTypeSlice []Read, fromDB bool) [][]Read {
	var currCell string
	var bin, count int
	var found bool

	whichBin := 'A'                   //for counting cells per bin
	whichBinMap := make(map[rune]int) //for counting cells per bin

	binnedCells := make([][]Read, s.BinCells) //make a slice of slice of string with the size of the user-specified number of bins
	prob := 1.0 / float64(s.BinCells)         // determine the probability that the cell belongs to a particular bin

	SortReadByCellBx(cellTypeSlice) //sort the slice of strings containing (cellBarcode \t construct) so that indentical cell barcodes line up next to one another
	for _, i := range cellTypeSlice {
		if i.Construct == "" {
			fmt.Println("distribute cells empty")
		}
		if i.Bx == currCell { // if the cell barcode is the same as the one the loop just saw, partition that cell-construct pair into the same bin as the one before
			binnedCells[bin] = append(binnedCells[bin], i)
			continue
		}
		bin = DetermineBin(s.BinCells, prob)           //function to determine which bin to put the new cell into.
		binnedCells[bin] = append(binnedCells[bin], i) //put cell into correct bin
		currCell = i.Bx                                //set the cell barcode that was just partitioned to current cell for the next iteration
		count, found = whichBinMap[whichBin+rune(bin)]
		if !found {
			whichBinMap[whichBin+rune(bin)] = 1
		} else {
			whichBinMap[whichBin+rune(bin)] = count + 1
		}
	}

	var outSlice []string
	if !fromDB {
		for i := range whichBinMap {
			outSlice = append(outSlice, fmt.Sprintf("Bin: %c\tcells: %d", i, whichBinMap[i]))
		}
		sort.Strings(outSlice)
		for _, i := range outSlice {
			fmt.Println(i)
		}
	}
	return binnedCells
}

// DetermineIdealBins will determine how many bins should be used for the pseudobulk -binCells option.
func DetermineIdealBins(s ScStarrSeqSettings, umiSlice []Read) int {
	var count, j int
	var found, stop bool
	var z, m string
	var binSlice []int
	var matrix [][]Read
	var l Read

	nc := fileio.Read(s.DetermineBins)

	for i := 0; i <= 10; i++ { //do this whole loop 10 times (essentially 10 replicates)
		stop = false
		for j = 1; j <= 100; j++ { //loop from bins 1 to 100 (j)
			s.BinCells = j                              //set the binCells parameter to j
			matrix = DistributeCells(s, umiSlice, true) //run distributeCells with that value of j
			for k := range matrix {                     //loop through bins in matrix
				ncMap := make(map[string]int) //create a map that has ncConstructName--counts
				for _, z = range nc {         //populate the map
					ncMap[z] = 0
				}
				for _, l = range matrix[k] { //loop through the constructs in a particular bin
					count, found = ncMap[l.Construct] //search the 2nd column to see if it is a negative control
					if found {                        //if it is, add a count to the map for that particular negative control
						ncMap[l.Construct] = count + 1
					}
				}
				for m = range ncMap {
					count, _ = ncMap[m]
					if count < 1 {
						binSlice = append(binSlice, j-1)
						stop = true
						break
					}
				}
				if stop == true {
					break
				}
			}
			if stop == true {
				break
			}
		}
	}
	return int(numbers.AverageInt(binSlice))
}

// BinnedPseudobulk is the psuedobulk function if the -binCells option is used. It adds an addition column to the dataframe corresponding to bin identity.
// It will handle writing of the pseudobulk maps to the provided EasyWriter input.
func BinnedPseudobulk(s ScStarrSeqSettings, inSlices [][]Read, out *fileio.EasyWriter) {
	var i string
	var count float64
	var found bool
	var j Read

	whichBin := 'A'
	fileio.WriteToFileHandle(out, "construct\tcounts\tbin")
	for _, bin := range inSlices {
		mp := make(map[string]float64)
		var toWrite []string
		for _, j = range bin {
			if j.Construct == "" {
				fmt.Println("empty")
			}
			count, found = mp[j.Construct]
			if !found {
				mp[j.Construct] = 1
			} else {
				mp[j.Construct] = count + 1
			}
		}

		if s.InputNormalize != "" {
			InputNormalize(mp, s.InputNormalize)
		}
		for i = range mp {
			toWrite = append(toWrite, fmt.Sprintf("%s\t%f\t%c", i, mp[i], whichBin))
		}
		sort.Strings(toWrite)
		for _, i = range toWrite {
			fileio.WriteToFileHandle(out, i)
		}
		whichBin++
	}
}
