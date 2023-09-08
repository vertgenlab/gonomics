package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"math/rand"
	"sort"
	"strings"
)

type ScStarrSeqSettings struct {
	InFile           string
	OutFile          string
	InputNormalize   string
	ByCell           string
	SamOut           string
	ScAnalysis       string
	BinCells         int
	UmiSat           string
	TransfectionNorm string
	Bed              string
	NcNorm           string
	DetermineBins    string
	InputSequencing  string
	SetSeed          int64
	CountMatrix      string
	NoOut            bool
}

// UMI is a custom struct that stores data about a UMI-representitive STARR-seq read. It stores the cell barcode, the cell type and the construct that the read maps to, all as strings
type UMI struct {
	Bx        string
	Cluster   string
	Construct string
}

// ClusterKey is a custom struct that stores a cell barcode and the cell type that that barcode belongs to, both as strings
type ClusterKey struct {
	Cluster string
	Bx      string
}

// SortUmiByCellBx takes in a slice of UMI and sorts the slice lexicographically according to the cell barcode
func SortUmiByCellBx(umiSlice []UMI) {
	sort.Slice(umiSlice, func(i, j int) bool {
		return umiSlice[i].Bx < umiSlice[j].Bx
	})
}

// ReadClusterKey takes in a string corresponding to a tab-delimined text file that has cell barcode in the first column and cell type in the second column and returns a slice of ClusterKey
func ReadClusterKey(inFile string) []ClusterKey {
	var clusterList []ClusterKey
	var clusterInfo ClusterKey
	var columns []string
	pairs := fileio.Read(inFile)
	for _, i := range pairs {
		columns = strings.Split(i, "\t")
		clusterInfo = ClusterKey{Bx: columns[0], Cluster: columns[1]}
		clusterList = append(clusterList, clusterInfo)
	}
	return clusterList
}

// ReadInputNormTable takes in a string corresponding to a tab-delimined text file that has constructName in the first column and input normalization factor
// in the second column and returns a slice of ClusterKey
func ReadInputNormTable(inFile string) []InputNormFactor {
	var inputNormList []InputNormFactor
	var inputNormInfo InputNormFactor
	var columns []string
	pairs := fileio.Read(inFile)
	for _, i := range pairs {
		columns = strings.Split(i, "\t")
		inputNormInfo = InputNormFactor{Construct: columns[0], NormFactor: parse.StringToFloat64(columns[1])}
		inputNormList = append(inputNormList, inputNormInfo)
	}
	return inputNormList
}

// UmiSaturation randomly subsets the whole bam file (10% to 100% of all reads) and calculates how many UMIs are in those subests. The output is a tab delimited text file.
func UmiSaturation(umiBxSlice []string, file string) {
	var perc, randNum float64
	var j string
	var count int

	out := fileio.EasyCreate(file)
	fileio.WriteToFileHandle(out, "totalReads\tumis")

	for i := 1; i <= 10; i++ {
		perc = float64(i) / 10
		mp := make(map[string]bool)
		count = 0
		for _, j = range umiBxSlice {
			randNum = rand.Float64()
			if randNum > perc {
				continue
			}
			count++
			mp[j] = true
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%d", count, len(mp)))
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func sortScCountMatrixByCluster(in []scStarrSeqMatrix) {
	sort.Slice(in, func(i, j int) bool {
		return in[i].cluster < in[j].cluster
	})
}

func sortScCountMatrixByConstruct(in []scStarrSeqMatrix) {
	sort.Slice(in, func(i, j int) bool {
		return in[i].construct < in[j].construct
	})
}
