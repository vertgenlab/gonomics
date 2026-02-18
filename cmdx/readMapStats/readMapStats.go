package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

type settings struct {
	selfparalogy string //flattened file
	peaks        string //peak set file
	alnfile      string //bam file of alignments
	outfile      string
}

func readMapStats(s settings) {
	var paraAns, peakAns []interval.Interval

	mp := make(map[string]map[uint8]int)
	mp["allSP"] = make(map[uint8]int)
	mp["sp_peaks"] = make(map[uint8]int)
	mp["noSP"] = make(map[uint8]int)
	mp["noSP_peaks"] = make(map[uint8]int)

	selfpara := bed.Read(s.selfparalogy) //read in the flattened self-paralogy file
	peakset := bed.Read(s.peaks)         //read in the peak file

	paraTree := interval.BedSliceToIntervalMap(selfpara)
	peakTree := interval.BedSliceToIntervalMap(peakset)

	ch, _ := sam.GoReadToChan(s.alnfile)
	for i := range ch {
		paraAns = interval.Query(paraTree, i, "any")
		peakAns = interval.Query(peakTree, i, "any")

		switch {
		case len(paraAns) > 0:
			mp["allSP"][i.MapQ]++
			if len(peakAns) > 0 {
				mp["sp_peaks"][i.MapQ]++
			}
		case len(paraAns) == 0:
			mp["noSP"][i.MapQ]++
			if len(peakAns) > 0 {
				mp["noSP_peaks"][i.MapQ]++
			}
		default:
			log.Fatalf("should get here")
		}

	}
	printMap(mp, s)

}

func printMap(mp map[string]map[uint8]int, s settings) {
	var j int
	var i uint8
	var ans []int
	var ansStr string

	out := fileio.EasyCreate(s.outfile)

	fileio.WriteToFileHandle(out, "mapq,allSP,SP_peaks,noSP,noSP_peaks")

	maps := []string{"allSP", "sp_peaks", "noSP", "noSP_peaks"}

	for i = 0; i < 61; i++ {
		ans = ans[:0]
		for j = range maps {
			ans = append(ans, mp[maps[j]][i])
		}
		ansStr = fileio.IntSliceToString(ans)
		fileio.WriteToFileHandle(out, fmt.Sprintf("%d,%s", i, ansStr))
	}
	exception.PanicOnErr(out.Close())
}

func main() {

	flag.Parse() // arguments are as follows: 0 -- aln file; 1 -- self paralogy file; 2 -- peak file; 3 -- outfile

	s := settings{
		selfparalogy: flag.Arg(1),
		peaks:        flag.Arg(2),
		alnfile:      flag.Arg(0),
		outfile:      flag.Arg(3),
	}
	readMapStats(s)
}

/*
func readMapStats(s settings) {
	var inSP []sam.Sam

	selfpara := bed.Read(s.selfparalogy) //read in the flattened self-paralogy file
	peakset := bed.Read(s.peaks)         //read in the peak file

	br, _ := sam.OpenBam(s.bamfile)
	bai := sam.ReadBai(s.bamfile + ".bai") //must have the .bai extension and be in the same dir

	for i := range selfpara {
		inSP = sam.SeekBamRegion(br, bai, selfpara[i].Chrom, uint32(selfpara[i].ChromStart), uint32(selfpara[i].ChromEnd))
	}


}

*/
