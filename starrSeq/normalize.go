package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

// InputNormFactor is a custom struct that scores input-normalization factors for STARR-seq constructs
type InputNormFactor struct {
	Construct  string
	NormFactor float64
}

// GfpNormalize iterates over a single-cell counts map and applies a GFP normalization factor.
func GfpNormalize(countsMap map[string]float64, normFactor float64) {
	var counts float64
	for i := range countsMap {
		counts, _ = countsMap[i]
		countsMap[i] = counts * normFactor
	}
}

// CalculateGfpNormFactor takes in a map of cellType-gfpReads and returns a map with cellType-gfpNormalizationValue. Normalization value is calculated by the percentage of GFP reads if all
// cell types were equal, and dividing it by the observed percentage of GFP reads in a given cell type
func CalculateGfpNormFactor(clusterGFP map[string]int) map[string]float64 {
	var totalCounts, clusterCounts int
	var actualPerc float64
	gfpNormMap := make(map[string]float64)
	idealPerc := 1.0 / float64(len(clusterGFP))
	for i := range clusterGFP {
		clusterCounts, _ = clusterGFP[i]
		totalCounts = clusterCounts + totalCounts
	}
	for i := range clusterGFP {
		actualPerc = float64(clusterGFP[i]) / float64(totalCounts)
		gfpNormMap[i] = idealPerc / actualPerc
	}
	return gfpNormMap
}

// GetNcAvg takes in the input normalized matrix of single cell data, the list of negative control sequences and an empty slice of slices and will return the average negative control reads for each
// cell type
func GetNcAvg(matrix []scStarrSeqMatrix, ncSlice []string, ncVals [][]float64) []float64 {
	var found bool
	var currCluster string
	var prevCluster string
	var ncAvg []float64

	clusterIdx := 0
	for p, i := range matrix {
		found = false
		currCluster = i.cluster
		if p > 0 && currCluster != prevCluster {
			clusterIdx++
		}
		for _, j := range ncSlice {
			if j == i.construct {
				found = true
				break
			}
		}
		if found {
			ncVals[clusterIdx] = append(ncVals[clusterIdx], i.counts)
		}
		prevCluster = currCluster
	}
	for i := range ncVals {
		ncAvg = append(ncAvg, numbers.AverageFloat64(ncVals[i]))
	}
	return ncAvg
}

// InputNormalize takes in the psuedobulk map and an input normalization table and normalizes all the raw count values in the map
func InputNormalize(mp map[string]float64, normalize string) {
	var total float64
	var found bool

	inputNormValues := ReadInputNormTable(normalize)

	if len(inputNormValues) < len(mp) {
		fmt.Println("The input normalization table has less constructs than were found in the input bam. 1 or more constructs won't be normalized. Please check your input files.")
	} else if len(inputNormValues) > len(mp) {
		fmt.Println("The input normalization table has more constructs than were found in the input bam. Constructs not found in the bam will have 0 counts in the output.")
	}
	for _, i := range inputNormValues { //iterate over the table with input normalization values and use those to edit the counts in the map
		total, found = mp[i.Construct]
		if found {
			mp[i.Construct] = total * i.NormFactor
		} else {
			mp[i.Construct] = 0.0 //user wants to normalize this construct meaning that it was in the library. however, no reads were recovered for that construct so we add it to the map and set the value to zero
		}
	}
}

// NormScToNegativeCtrls is the initial function to normalize single-cell data to negative control reads in the same cell type. It takes in a matrix that is the input-normalized, the settings struct,
// the number of cell types and returns the negative control normalized read counts
func NormScToNegativeCtrls(matrix []scStarrSeqMatrix, ncNorm string, numCellTypes int) []scStarrSeqMatrix {
	var ncSlice []string
	var prevCluster, currCluster string

	ncVals := make([][]float64, numCellTypes)
	nc := fileio.Read(ncNorm)
	for _, i := range nc {
		ncSlice = append(ncSlice, i)
	}
	ncAvg := GetNcAvg(matrix, ncSlice, ncVals)

	clusterIdx := 0
	for j, i := range matrix {
		currCluster = i.cluster
		if j > 0 && currCluster != prevCluster {
			clusterIdx++
		}
		i.counts = i.counts / ncAvg[clusterIdx]
		prevCluster = currCluster
	}
	return matrix
}

// ParseGfpBam is similar to the parseBam function but handles bams containing GFP reads. This function also takes in a map of cellBx-cellType and an empty map containing each cellType
// It returns a map that contains [cellType] = gfpReads
func ParseGfpBam(s ScStarrSeqSettings, cellTypeMap map[string]string, clusterGFP map[string]int) ([]Read, map[string]int) {
	var bit int32
	var cluster, cellBxFormat string
	var gfpUmis []Read
	var gfpStats []string
	var found bool
	var multipleGems bool = false
	var gemNumber int = 1

	//detect multiple bam files
	files := strings.Split(s.TransfectionNorm, ",")
	if len(files) > 1 {
		multipleGems = true
		if s.ScAnalysis != "" {
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -transfectionNorm. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		}
	}

	for _, bam := range files {
		inChan, _ := sam.GoReadToChan(bam)
		for i := range inChan {
			num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
			bit = num.(int32)
			if bit&8 == 8 {
				cellBx, _, _ := sam.QueryTag(i, "CB")
				if multipleGems {
					cellBxFormat = fmt.Sprintf("%s_%d", cellBx.(string), gemNumber)
				} else {
					cellBxFormat = cellBx.(string)
				}
				cluster, found = cellTypeMap[cellBxFormat]
				switch {
				case found:
					gfpUmis = append(gfpUmis, Read{Bx: cellBxFormat, Cluster: cluster, Construct: "GFP"})
					clusterGFP[cluster] += 1
				case !found:
					gfpUmis = append(gfpUmis, Read{Bx: cellBxFormat, Cluster: "undefined", Construct: "GFP"})
				}
			}
		}
		gemNumber++
	}

	//print out some GFP read counts / cluster
	var gfpReads int
	for i := range clusterGFP {
		gfpReads, _ = clusterGFP[i]
		gfpStats = append(gfpStats, fmt.Sprintf("found %d GFP reads in %s", gfpReads, i))
	}
	for _, i := range gfpStats {
		fmt.Println(i)
	}
	return gfpUmis, clusterGFP
}

// MakeInputNormMap takes in a file that has tab-delimited "construct" and "input normalization factor" pairs and returns a map of those pairs
func MakeInputNormMap(inputNormTable string) map[string]float64 {
	mp := make(map[string]float64)
	file := ReadInputNormTable(inputNormTable)
	for _, i := range file {
		mp[i.Construct] = i.NormFactor
	}
	return mp
}
