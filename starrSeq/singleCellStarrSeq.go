package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
)

type scStarrSeqMatrix struct {
	cluster   string
	construct string
	counts    float64
}

// ClusterKey is a custom struct that stores a cell barcode and the cell type that that barcode belongs to, both as strings
type ClusterKey struct {
	Cluster string
	Bx      string
}

// SingleCellAnalysis will create a count for each construct in each cell type. This function takes in starrSeqSetting, a readSlice and some cell type info. It will create
// an output file that has readCounts for each construct in each cell type present. It is compatible with input normalization and either transfection reporter transfection
// cluster normalization or negative control cluster normalization
func SingleCellAnalysis(s ScStarrSeqSettings, readSlice []Read, allConstructs []string, allCellTypes []string, cellTypeMap map[string]string, out *fileio.EasyWriter) {
	var cellType, j, m, toWrite string
	var outMatrix []scStarrSeqMatrix
	var outLine scStarrSeqMatrix
	var gfpNormFactorMap map[string]float64

	allCellTypesMap := make(map[string]int) // [cellType]placeholderInt
	ncMap := make(map[string]float64)       // [ncConstruct]Counts

	for _, i := range cellTypeMap {
		allCellTypesMap[i] = 0
	}

	if s.TransfectionNorm != "" {
		_, gfpClusterMap := ParseGfpBam(s, cellTypeMap, allCellTypesMap)
		gfpNormFactorMap = CalculateGfpNormFactor(gfpClusterMap)
	}

	if s.NcNorm != "" {
		ncNames := fileio.Read(s.NcNorm)
		for _, i := range ncNames {
			ncMap[i] = 0
		}
	}

	for cellType = range allCellTypesMap { //loop over the list of cell types
		var a int = 0
		singleCellTypeMap := make(map[string]float64) // make a map which will hold the construct counts for the current cell type ([construct]counts)
		for _, j = range allConstructs {              // loop over all construct that were found in the bam file to make every construct start at 0 reads
			singleCellTypeMap[j] = 0
		}
		for _, i := range readSlice { //loop over the slice created in cellrangerBam() that that contains: cellBarcode \t construct
			if i.Cluster == cellType { // if that cell is found and the cell type in the map matches the cell type we are looping through, search for the construct corresponding with that cellBarcode in the singleCellCount map and add one to the value
				singleCellTypeMap[i.Construct] += 1
				a++
			}
		}
		if s.InputNormalize != "" { // normalize the singleCellCount map if needed
			InputNormalize(singleCellTypeMap, s.InputNormalize)
		}
		if s.TransfectionNorm != "" {
			gfpNormalizationFactor := gfpNormFactorMap[cellType]
			GfpNormalize(singleCellTypeMap, gfpNormalizationFactor)
		}
		for m = range singleCellTypeMap { // write out the results from that cell type
			outLine = scStarrSeqMatrix{counts: singleCellTypeMap[m], cluster: cellType, construct: m}
			outMatrix = append(outMatrix, outLine)
		}
		fmt.Println(fmt.Sprintf("Found %d raw counts in the cell type: %s", a, cellType))
	}

	sortScCountMatrixByConstruct(outMatrix)
	sortScCountMatrixByCluster(outMatrix) //sort the slice before writing

	if s.NcNorm != "" {
		outMatrix = NormScToNegativeCtrls(outMatrix, s.NcNorm, len(allCellTypes))
	}
	if s.NcNorm != "" {
		fileio.WriteToFileHandle(out, "cellCluster\tconstruct\tcounts/ncCounts")
	} else {
		fileio.WriteToFileHandle(out, "cellCluster\tconstruct\tcounts")
	}

	for _, i := range outMatrix { //once all cell types have been looped through, write out the total count matrix
		toWrite = fmt.Sprintf("%s\t%s\t%f", i.cluster, i.construct, i.counts)
		fileio.WriteToFileHandle(out, toWrite)
	}
}
