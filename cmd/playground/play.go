package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
)

type test struct {
	i int
	s string
	b bool
	p *test
}

func main() {
	/*chroms := []string{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
		"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
		"chr22", "chrX", "chrY"}

	for i := range chroms {
		fmt.Printf("lastz ../hs1.byChrom/%s.hs1.fa --self --output=withMask/axtTest/selfAlign/%s_%s.hs1.lav --allocate:traceback=1G --scores=human_chimp_v2.mat O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0\n", chroms[i], chroms[i], chroms[i])
		//for j := i + 1; j < len(chroms); j++ {
		//	fmt.Printf("lastz ../hs1.byChrom/%s.hs1.fa ../hs1.byChrom/%s.hs1.fa --format=axt --output=withMask/axtTest/%s_%s.hs1.axt --allocate:traceback=1G --scores=human_chimp_v2.mat O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0\n", chroms[i], chroms[j], chroms[i], chroms[j])
		//}
	}


		var cols []string
		txt := fileio.Read("/Users/sethweaver/Downloads/hsSD/netSdDiscovery/pt.chromNames.txt")
		for i := range txt {
			cols = strings.Split(txt[i], "\t")
			fmt.Printf("sed -i 's/%s/%s/g hs1.chainSynGCA_028858775.2.bed \n", cols[0], cols[1])
		}


	*/

	/*
			var fd, fd2 familyData
			var cols, extraInfo []string
			var rs refSeqData
			var found, found2 bool
			var err error

			ncbiRefSeq := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/hs1.ncbiRefSeq.bed")
			hgnc := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/HGNC/HGNC.geneNamesWithFamily.txt")
			outGood := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/HGNC/hs1.ncbiRefSeq.familyData.bed")
			outPartial := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/HGNC/hs1.ncbiRefSeq.ambigousFamData.txt")
			outNotFound := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/HGNC/hs1.ncbiRefSeq.NoFamData.bed")
			nameMap, refSeqMap := make(map[string]familyData), make(map[string]familyData)

			for i := range hgnc {
				cols = strings.Split(hgnc[i], "\t")
				fd = familyData{
					HGNC_ID:   cols[0],
					geneName:  cols[1],
					RefSeqID:  cols[2],
					GeneGroup: cols[4],
					ensemblID: cols[3],
				}
				nameMap[fd.geneName] = fd
				refSeqMap[fd.RefSeqID] = fd
			}

			for i := range ncbiRefSeq {
				cols = strings.Split(ncbiRefSeq[i], "\t")
				extraInfo = strings.Split(cols[4], " ")
				rs = refSeqData{
					chrom:      cols[0],
					chromStart: cols[1],
					chromEnd:   cols[2],
					strand:     cols[3],
					geneName:   extraInfo[1][1 : len(extraInfo[1])-2],
					refSeqId:   extraInfo[3][1 : len(extraInfo[3])-4],
				}
				fd, found = nameMap[rs.geneName]
				fd2, found2 = refSeqMap[rs.refSeqId]
				if found || found2 {
					if (fd.HGNC_ID == fd2.HGNC_ID) && fd.HGNC_ID != "" {
						rs.EnsemblID = fd.ensemblID
						rs.familyID = fd.GeneGroup
						rs.HGNC_ID = fd.HGNC_ID
						writeGoodFind(outGood, rs)
					} else {
						writePartial(outPartial, ncbiRefSeq[i], fd, fd2)
					}
				} else {
					fileio.WriteToFileHandle(outNotFound, ncbiRefSeq[i])
				}
			}
			err = outGood.Close()
			exception.PanicOnErr(err)
			err = outPartial.Close()
			exception.PanicOnErr(err)
			err = outNotFound.Close()
			exception.PanicOnErr(err)
		}

		func writeGoodFind(outGood *fileio.EasyWriter, rs refSeqData) {
			fileio.WriteToFileHandle(outGood, fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", rs.chrom, rs.chromStart, rs.chromEnd, rs.strand, rs.geneName, rs.refSeqId, rs.EnsemblID, rs.HGNC_ID, rs.familyID))
		}

		func writePartial(outPartial *fileio.EasyWriter, bedLine string, fd, fd2 familyData) {
			fileio.WriteToFileHandle(outPartial, strings.Join([]string{bedLine, fdToString(fd), fdToString(fd2)}, "\t"))
		}

		func fdToString(fd familyData) string {
			return fmt.Sprintf("%s\t%s\t%s\t%s\t%s", fd.HGNC_ID, fd.geneName, fd.RefSeqID, fd.ensemblID, fd.GeneGroup)

	*/

	/*
		var cols []string
		var rs refSeqData

		txt := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/HGNC/hs1.ncbiRefSeq.ambigousButSaved.bed")

		out := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/HGNC/hs1.ncbiRefSeq.ambigousButSaved.bed.tmp")

		for i := range txt {
			cols = strings.Split(txt[i], "\t")
			rs = refSeqData{
				chrom:      cols[0],
				chromStart: cols[1],
				chromEnd:   cols[2],
				strand:     cols[3],
				geneName:   cols[4],
				HGNC_ID:    cols[5],
			}
			if strings.Contains(cols[6], "_") {
				rs.refSeqId = cols[6]
				if strings.Contains(cols[7], "ENSG") {
					rs.EnsemblID = cols[7]
					rs.familyID = cols[8]
				} else {
					rs.EnsemblID = ""
					rs.familyID = cols[7]
				}
			} else {
				rs.refSeqId = ""
				if strings.Contains(cols[6], "ENSG") {
					rs.EnsemblID = cols[6]
					rs.familyID = cols[7]
				} else {
					rs.EnsemblID = ""
					rs.familyID = cols[6]
				}
			}
			fileio.WriteToFileHandle(out, rsToString(rs))
		}
		exception.PanicOnErr(out.Close())

	*/

	/*
		var cols []string
		var gtexData string
		var found1, found2 bool
		var c int

		gtexMapES, gtexMapGN := make(map[string]string), make(map[string]string)
		gtex := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/gtexHg38.txt")
		out := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/genesInHs1SegDup.withFamData.overlapThresh1.withGTEX.bed")

		for i := range gtex {
			cols = strings.Split(gtex[i], "\t")
			gtexMapES[strings.Split(cols[6], ".")[0]] = gtex[i]
			gtexMapGN[cols[3]] = gtex[i]
		}
		//fmt.Println(gtexMap)

		genesInSD := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/genesInHs1SegDup.withFamData.overlapThresh1.bed")

		for i := range genesInSD {
			cols = strings.Split(genesInSD[i], "\t")
			gtexData, found1 = gtexMapES[cols[6]]
			if found1 {
				fileio.WriteToFileHandle(out, strings.Join([]string{genesInSD[i], gtexData}, "\t"))
				continue
			}
			gtexData, found2 = gtexMapGN[cols[4]]
			if found2 {
				fileio.WriteToFileHandle(out, strings.Join([]string{genesInSD[i], gtexData}, "\t"))
				continue
			}
			c++
		}
		fmt.Println(c)
		exception.PanicOnErr(out.Close())

	*/
	var cols, fams, slc []string
	var gex, gexByTissue [][]float64
	var j int
	var sd float64
	var stddev []float64

	data := fileio.Read("/Users/sethweaver/Downloads/hsSD/genes/genesInHs1SegDup.withFamData.overlapThresh1.withGTEX.bed")
	out := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/genes/gtexStdDevSdGeneFams.txt")

	mp := make(map[int][]string)

	for i := range data {
		cols = strings.Split(data[i], "\t")
		fams = strings.Split(cols[8], "|")
		if fams[0] == "" {
			continue
		}
		for j = range fams {
			mp[parse.StringToInt(fams[j])] = append(mp[parse.StringToInt(fams[j])], data[i])
		}
	}

	fileio.WriteToFileHandle(out, "familyID\tnumGenes\tstddevPerTissue")
	for i := range mp {
		stddev = []float64{}
		gex = [][]float64{}
		slc = mp[i]
		if len(slc) == 1 {
			continue
		}
		for j = range slc {
			cols = strings.Split(slc[j], "\t")
			gex = append(gex, fileio.StringToFloatSlice(cols[18]))
		}
		gexByTissue = geneToTissue(gex)
		for j = range gexByTissue {
			sd = numbers.StandardDeviationFloat64(gexByTissue[j])
			stddev = append(stddev, sd)
		}
		fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%d\t%s", i, len(slc), fileio.FloatSliceToString(stddev)))
	}
	exception.PanicOnErr(out.Close())
}

func geneToTissue(gex [][]float64) [][]float64 {
	var gexByTissue [][]float64
	for i := 0; i < len(gex[0]); i++ {
		gexByTissue = append(gexByTissue, make([]float64, len(gex)))
		for j := 0; j < len(gex); j++ {
			gexByTissue[i][j] = gex[j][i]
		}
	}
	return gexByTissue
}

func rsToString(rs refSeqData) string {
	return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", rs.chrom, rs.chromStart, rs.chromEnd, rs.strand, rs.geneName, rs.refSeqId, rs.EnsemblID, rs.HGNC_ID, rs.familyID)
}

type familyData struct {
	HGNC_ID   string
	geneName  string
	RefSeqID  string
	GeneGroup string
	ensemblID string
}

type refSeqData struct {
	chrom      string
	chromStart string
	chromEnd   string
	strand     string
	geneName   string
	refSeqId   string
	EnsemblID  string
	familyID   string
	HGNC_ID    string
}
