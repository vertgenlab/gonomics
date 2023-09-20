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
	InFile               string
	OutFile              string
	InputNormalize       string
	ByCell               string
	SamOut               string
	ScAnalysis           string
	BinCells             int
	UmiSat               string
	TransfectionNorm     string
	Bed                  string
	NcNorm               string
	DetermineBins        string
	InputSequencing      string
	SetSeed              int64
	CountMatrix          string
	NoOut                bool
	AltMapping           string
	CountMatrixCellTypes string
}

type ReadSliceAnalysisSettings struct {
	FuncSettings  ScStarrSeqSettings
	UmiBxSlice    []string
	AllCellTypes  []string
	AllConstructs []string
	ReadSlice     []Read
	CellTypeMap   map[string]string
}

// Read is a custom struct that stores data about a UMI-representitive STARR-seq read. It stores the cell barcode, the cell type and the construct that the read maps to, all as strings
type Read struct {
	Bx        string
	Cluster   string
	Construct string
	UMI       string
}

// ClusterKey is a custom struct that stores a cell barcode and the cell type that that barcode belongs to, both as strings
type ClusterKey struct {
	Cluster string
	Bx      string
}

// SortUmiByCellBx takes in a slice of Read and sorts the slice lexicographically according to the cell barcode
func SortUmiByCellBx(umiSlice []Read) {
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

// UmiSaturation randomly subsets the whole bam file (10% to 100% of all reads) and calculates how many Reads are in those subests. The output is a tab delimited text file.
func UmiSaturation(s ScStarrSeqSettings, umiBxSlice []string) {
	var perc, randNum float64
	var j string
	var count int

	out := fileio.EasyCreate(s.UmiSat)
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

// WritePseudobulkMap writes out a pseudobulk map to an io.writer
func WritePseudobulkMap(mp map[string]float64, writer *fileio.EasyWriter) {
	var total float64
	var write string
	var writeSlice []string
	for i := range mp {
		total, _ = mp[i]
		write = fmt.Sprintf("%s\t%f", i, total)
		writeSlice = append(writeSlice, write)
	}
	sort.Strings(writeSlice)
	fileio.WriteToFileHandle(writer, "construct\tcounts")
	for _, i := range writeSlice {
		fileio.WriteToFileHandle(writer, i)
	}
}

// ReadSliceAnalysis takes in a slice of Read and the settings struct and distributes the Read slice to the specified analysis functions
func ReadSliceAnalysis(r ReadSliceAnalysisSettings) {
	var out *fileio.EasyWriter
	//all non-defult output options
	if r.FuncSettings.ByCell != "" {
		ByCell(r.FuncSettings, r.ReadSlice)
	}
	if r.FuncSettings.UmiSat != "" {
		UmiSaturation(r.FuncSettings, r.UmiBxSlice)
	}
	if r.FuncSettings.CountMatrix != "" {
		MakeCountMatrix(r.FuncSettings, r.ReadSlice, r.AllConstructs)
	}
	// default output options
	if !r.FuncSettings.NoOut {
		out = fileio.EasyCreate(r.FuncSettings.OutFile)
		if r.FuncSettings.ScAnalysis != "" {
			SingleCellAnalysis(r.FuncSettings, r.ReadSlice, r.AllConstructs, r.AllConstructs, r.CellTypeMap, out)
		} else if r.FuncSettings.DetermineBins != "" {
			r.FuncSettings.BinCells = DetermineIdealBins(r.FuncSettings, r.ReadSlice)
			fmt.Printf("Using %d bins\n", r.FuncSettings.BinCells)
		}
		if r.FuncSettings.BinCells > 1 {
			binnedCells := DistributeCells(r.FuncSettings, r.ReadSlice, false)
			BinnedPseudobulk(r.FuncSettings, binnedCells, out)
		}
		if r.FuncSettings.ScAnalysis == "" && r.FuncSettings.BinCells <= 1 {
			//default
			pbMap := ReadSliceToPseudobulk(r.FuncSettings, r.ReadSlice)
			WritePseudobulkMap(pbMap, out)
		}
		err := out.Close()
		exception.PanicOnErr(err)
	}
}

// ReadSliceToPseudobulk takes a settings struct slice of Read and returns a map of construct -- read counts. If an input normalization
// table is present in the settings struct. The returned map will be input normalized
func ReadSliceToPseudobulk(s ScStarrSeqSettings, readSlice []Read) map[string]float64 {
	constructMap := make(map[string]float64)

	for _, i := range readSlice {
		constructMap[i.Construct] += 1
	}
	if s.InputNormalize != "" {
		InputNormalize(constructMap, s.InputNormalize)
	}
	return constructMap
}

func ByCell(s ScStarrSeqSettings, readSlice []Read) {
	outByCell := fileio.EasyCreate(s.ByCell)
	for _, i := range readSlice {
		fileio.WriteToFileHandle(outByCell, fmt.Sprintf("%s\t%s", i.Bx, i.Construct))
	}
	err := outByCell.Close()
	exception.PanicOnErr(err)
}
