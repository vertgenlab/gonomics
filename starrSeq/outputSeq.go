package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"sort"
	"strings"
)

// parseBam takes in a cellranger count bam file and pulls out reads that are representative of the UMI and also returns the construct associated with the UMI.
func ParseCellrangerBam(s OutputSeqSettings) {
	var k int = 0
	var constructName, umiBx, cellType, cellFormat string
	var found bool
	var bit int32
	var multipleGems = false
	var gemNumber int = 1
	var allConstructs, umiBxSlice, allCellTypes []string
	var read Read
	var readSlice []Read
	var outSam *fileio.EasyWriter
	var ck []ClusterKey

	//detect multiple bam files
	files := strings.Split(s.InFile, ",")
	if len(files) > 1 {
		multipleGems = true
		if s.ScAnalysis != "" {
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -gfpNorm. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		}
	}

	cellTypeMap := make(map[string]string)
	allCellTypesMap := make(map[string]int)

	if s.ScAnalysis != "" || s.CountMatrixCellTypes != "" {
		if s.ScAnalysis != "" {
			ck = ReadClusterKey(s.ScAnalysis)
		} else {
			ck = ReadClusterKey(s.CountMatrixCellTypes)
		}
		for _, i := range ck {
			cellTypeMap[i.Bx] = i.Cluster
			allCellTypesMap[i.Cluster] = 0
		}
		for i := range allCellTypesMap {
			allCellTypes = append(allCellTypes, i)
		}
		sort.Strings(allCellTypes)
	}

	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	if s.Bed != "" {
		bedEntries := bed.Read(s.Bed)
		for _, i := range bedEntries {
			tree = append(tree, i)
		}
		selectTree = interval.BuildTree(tree)
	}

	for _, bam := range files {
		ch, head := sam.GoReadToChan(bam)
		if s.SamOut != "" {
			if gemNumber == 1 {
				outSam = fileio.EasyCreate(s.SamOut)
				sam.WriteHeaderToFileHandle(outSam, head)
			}
		}
		for i := range ch { //iterate of over the chanel of sam.Sam
			if s.UmiSat != "" {
				readUmi, foundUMI, _ := sam.QueryTag(i, "UB")
				cb, foundCb, _ := sam.QueryTag(i, "CB")

				if foundCb && foundUMI {
					umiBx = fmt.Sprintf("%s_%s", readUmi, cb)
					umiBxSlice = append(umiBxSlice, umiBx)
				}
			}
			num, _, _ := sam.QueryTag(i, "xf") //xf: extra flags (cellranger flags)
			bit = num.(int32)
			if bit&8 == 8 { // bit 8 is the flag for a UMI that was used in final count. I call these "valid" UMIs.
				k++
				if s.Bed != "" {
					overlappedBedEntry := interval.Query(selectTree, i, "any")
					if len(overlappedBedEntry) != 1 {
						continue
					}
					constructName = overlappedBedEntry[0].(bed.Bed).Name
				} else {
					construct, _, _ := sam.QueryTag(i, "GX") // get the construct associated with valid UMI
					constructName = construct.(string)
				}
				//rand.Seed(s.setSeed)
				cell, _, _ := sam.QueryTag(i, "CB") // get the cell barcode associated with valid UMI
				if multipleGems {
					cellFormat = fmt.Sprintf("%s_%d", cell.(string), gemNumber)
				} else {
					cellFormat = cell.(string)
				}
				readUmi, _, _ := sam.QueryTag(i, "UB") // get the UMI associated with the valid UMI
				cellType, found = cellTypeMap[cellFormat]
				switch {
				case (s.ScAnalysis != "" || s.CountMatrixCellTypes != "") && found: //we want cell type data and its found
					read = Read{Bx: cellFormat, Cluster: cellType, Construct: constructName, UMI: readUmi.(string)}
				case (s.ScAnalysis != "" || s.CountMatrixCellTypes != "") && !found: //we want cell type data, but we didn't find the cell type in the map
					read = Read{Bx: cellFormat, Cluster: "undefined", Construct: constructName, UMI: readUmi.(string)}
				case s.ScAnalysis == "" && s.CountMatrixCellTypes == "": // we don't care about cell type analysis. cluster is nil
					read = Read{Bx: cellFormat, Construct: constructName, UMI: readUmi.(string)}
				}
				allConstructs = append(allConstructs, read.Construct) //list of all constructs found in the bam
				readSlice = append(readSlice, read)
				if s.SamOut != "" { //write output as sam
					sam.WriteToFileHandle(outSam, i)
				}
			}
		}
		gemNumber++
	}
	fmt.Println("Found this many valid UMIs: ", k)
	r := ReadSliceAnalysisSettings{ReadSlice: readSlice, FuncSettings: s, AllCellTypes: allCellTypes, AllConstructs: allConstructs, CellTypeMap: cellTypeMap, UmiBxSlice: umiBxSlice}
	ReadSliceAnalysis(r)
}
