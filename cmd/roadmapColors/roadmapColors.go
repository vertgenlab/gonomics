package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

func makeMap() map[string]string {
	mp := make(map[string]string)
	mp["1_TssA"] = "120,188,97"
	mp["2_TssAFlnk"] = "44,85,48"
	mp["3_TxFlnk"] = "115,158,130"
	mp["4_Tx"] = "122,229,130"
	mp["5_TxWk"] = "122,229,130"
	mp["6_EnhG"] = "0,165,207"
	mp["7_Enh"] = "0,78,100"
	mp["8_ZNF/Rpts"] = "255,49,46"
	mp["9_Het"] = "81,80,82"
	mp["10_TssBiv"] = "255,200,87"
	mp["11_BivFlnk"] = "175,155,70"
	mp["12_EnhBiv"] = "243,167,18"
	mp["13_ReprPC"] = "219,43,57"
	mp["14_ReprPCWk"] = "219,43,57"
	mp["15_Quies"] = "188,171,174"
	return mp
}

func write(b bed.Bed, o *fileio.EasyWriter, color string) {
	fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d\t%d\t%s\t0\t.\t%d\t%d\t%s", b.Chrom, b.ChromStart, b.ChromEnd, b.Name, b.ChromStart, b.ChromEnd, color))
}

func main() {
	flag.Parse()
	filename := flag.Arg(0)
	prefix := filename[0 : len(filename)-4]
	o := fileio.EasyCreate(fmt.Sprintf("%s.colors.bed", prefix))
	fileio.WriteToFileHandle(o, fmt.Sprintf("track name=\"%s\" itemRgb=\"on\"", prefix))
	ch := bed.GoReadToChan(filename)
	mp := makeMap()
	for i := range ch {
		write(i, o, mp[i.Name])
	}
	err := o.Close()
	exception.PanicOnErr(err)
}
