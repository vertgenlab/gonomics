package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	sort2 "github.com/vertgenlab/gonomics/sort"
	"strings"
)

type line struct {
	bx      string
	cluster string
	vals    []float64
}

func MakeCountMatrix(s ScStarrSeqSettings, umiSlice []UMI, allConstructs []string) {
	var header []string
	var cell string
	var norm bool = false
	var normMap map[string]float64
	var matrixSlice []line
	var idx int
	var l line

	out := fileio.EasyCreate(s.CountMatrix)

	cellTypeMap := make(map[string]string) //bx--cellType

	header = append(header, "cellBx")
	if s.ScAnalysis != "" {
		header = append(header, "cluster")
	}

	uniqueConstructs := sort2.Unique(allConstructs)
	idxMap := make(map[string]int) //construct name and what place in the slice it goes
	for i, j := range uniqueConstructs {
		idxMap[j] = i
		header = append(header, j)
	}
	if s.TransfectionNorm != "" {
		header = append(header, "GFP")
	}
	headSmush := strings.Join(header, "\t")

	idxNormMap := make(map[int]float64)
	if s.InputNormalize != "" {
		norm = true
		normMap = MakeInputNormMap(s.InputNormalize)
		for i, j := range uniqueConstructs {
			idxNormMap[i] = normMap[j]
		}
	}

	fileio.WriteToFileHandle(out, headSmush)
	SortUmiByCellBx(umiSlice)
	cell = umiSlice[0].Bx

	l.vals = make([]float64, len(uniqueConstructs))
	var firstTime bool = true
	for i := range umiSlice {
		if cell == umiSlice[i].Bx {
			//add to current line
			if firstTime {
				firstTime = false
				l.bx = umiSlice[i].Bx
				l.cluster = umiSlice[i].Cluster
			}
			idx = idxMap[umiSlice[i].Construct]
			l.vals[idx] = l.vals[idx] + 1
		} else {
			//move on
			if norm {
				l = normLine(l, idxNormMap)
			}
			if s.TransfectionNorm == "" {
				writeLineToFileHandle(l, out, s)
			} else {
				matrixSlice = append(matrixSlice, l)
			}
			l = line{
				bx:      umiSlice[i].Bx,
				cluster: umiSlice[i].Cluster,
				vals:    nil,
			}
			l.vals = make([]float64, len(uniqueConstructs))
			cell = umiSlice[i].Bx
			idx = idxMap[umiSlice[i].Construct]
			l.vals[idx] = l.vals[idx] + 1
		}
	}
	//deal with last line
	if norm {
		l = normLine(l, idxNormMap)
	}
	if s.TransfectionNorm == "" {
		writeLineToFileHandle(l, out, s)
	} else {
		matrixSlice = append(matrixSlice, l)
	}

	if s.TransfectionNorm != "" {
		var found bool
		//add GFP to matrix
		placeHolderMap := make(map[string]int)
		gfpUmiSlice, _ := ParseGfpBam(s.TransfectionNorm, cellTypeMap, placeHolderMap)
		cellGfpMap := make(map[string]float64)
		for _, i := range gfpUmiSlice {
			_, found = cellGfpMap[i.Bx]
			if !found {
				cellGfpMap[i.Bx] = 1
			} else {
				cellGfpMap[i.Bx] += 1
			}
		}
		for _, i := range matrixSlice {
			_, found = cellGfpMap[i.bx]
			if !found {
				i.vals = append(i.vals, 0)
			} else {
				i.vals = append(i.vals, cellGfpMap[i.bx])
			}
			writeLineToFileHandle(i, out, s)
		}
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func writeLineToFileHandle(l line, out *fileio.EasyWriter, s ScStarrSeqSettings) {
	var floatSlice []string
	for _, i := range l.vals {
		floatSlice = append(floatSlice, fmt.Sprintf("%f", i))
	}
	valSmush := strings.Join(floatSlice, "\t")
	if s.ScAnalysis != "" {
		toWrite := fmt.Sprintf("%s\t%s\t%s", l.bx, l.cluster, valSmush)
		fileio.WriteToFileHandle(out, toWrite)
	} else {
		toWrite := fmt.Sprintf("%s\t%s", l.bx, valSmush)
		fileio.WriteToFileHandle(out, toWrite)
	}
}

func normLine(l line, idxNormMap map[int]float64) line {
	for i := range l.vals {
		l.vals[i] = idxNormMap[i] * l.vals[i]
	}
	return l
}
