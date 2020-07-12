package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func editFasta(inputFa string, outputFa string, bedfile string, reverse bool) {
	fa := fasta.Read(inputFa)
	regions := make(chan *bed.Bed)
	go bed.ReadToChan(bedfile, regions)
	chromHash := make(map[string][]*bed.Bed)
	for i := range regions {
		chromHash[i.Chrom] = append(chromHash[i.Chrom], i)
	}
	writer := fileio.EasyCreate(outputFa)
	defer writer.Close()
	var curr []*bed.Bed
	for _, chr := range fa {
		curr = chromHash[chr.Name]
		if len(curr) > 0 {
			if reverse {
				curr = bed.InvertRegions(chromHash[chr.Name], len(chr.Seq))
			}
			fasta.WriteFasta(writer, bed.EditFastaByRegion(chr, curr), 50)
		} else if len(curr) == 0 { //if chrom does not contain any beds, write fasta as is
			//log.Printf("Empty chrom: %s len=%d\n", chr.Name, len(chr.Seq))
			fasta.WriteFasta(writer, &fasta.Fasta{Name: chr.Name, Seq: chr.Seq}, 50)
			log.Printf("wrote chr len=%d\n", len(chr.Seq))
		}
	}
}

func usage() {
	fmt.Print(
		"editFasta - edit fasta sequence using a bed region file\n\n" +
			"Usage:\n" +
			"  editFasta [options] in.fa out.file\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	var expectedNumArgs int = 2

	var bedfile *string = flag.String("bed", "", "use `.bed` to edit fasta to desired genomic ranges")
	var invert = flag.Bool("invert", false, "edit reverse genomic ranges provided in bed file")
	flag.Parse()

	var inFile string = flag.Arg(0)
	var outFile string = flag.Arg(1)

	if expectedNumArgs != len(flag.Args()) {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	} else if *bedfile == "" {
		log.Fatalf("Error: must provide bed file to edit fasta....\n")
	} else {
		editFasta(inFile, outFile, *bedfile, *invert)
	}
}

/*
//TODO: Belongs in the fasta package, but bed and fasta becomes an illegal import cycle
func chrSplitByNs(chr *fasta.Fasta) []*fasta.Fasta {
	unGapped := bed.UngappedRegionsFromFa(chr)
	var answer []*fasta.Fasta = make([]*fasta.Fasta, len(unGapped))
	for i := 0; i < len(unGapped); i++ {
		answer[i] = &fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d", unGapped[i].Chrom, unGapped[i].ChromStart, unGapped[i].ChromEnd), Seq: chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]}
	}
	return answer
}
func unitfyNsAll(filename string, outFile string) {
	reader := fileio.EasyOpen(filename)
	writer := fileio.EasyCreate(outFile)
	defer reader.Close()
	defer writer.Close()

	//var ans *fasta.Fasta
	for fa, done := fasta.NextFasta(reader); !done; fa, done = fasta.NextFasta(reader) {
		fasta.WriteFasta(writer, unitfyNs(fa), 50)

		//ans = append(ans, unitfyNs(fa[i]))
	}

}
func unitfyNs(chr *fasta.Fasta) *fasta.Fasta {
	var ans *fasta.Fasta = &fasta.Fasta{Name: chr.Name}

	unGapped := bed.UngappedRegionsFromFa(chr)
	if len(unGapped) == 1 {
		log.Printf("Name=%s\n", chr.Name)
		ans.Seq = append(ans.Seq, chr.Seq...)
		return chr
	}
	for i := 0; i < len(unGapped)-1; i++ {
		//if len(chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]) < 100 {
			//log.Printf("Error: seqeunce too short...\nlen=%d, %s\n", len(chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]), dna.BasesToString(chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]))
		//	ans.Seq = append(ans.Seq, dna.CreateAllNs(100)...)
			//return nil
		//} else {
			ans.Seq = append(ans.Seq, chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]...)
			ans.Seq = append(ans.Seq, dna.CreateAllNs(100)...)
		//}

	}
	if len(unGapped) > 0 && len(chr.Seq[unGapped[len(unGapped)-1].ChromStart:unGapped[len(unGapped)-1].ChromEnd]) < 100 {
		ans.Seq = append(ans.Seq, chr.Seq[unGapped[len(unGapped)-1].ChromStart:unGapped[len(unGapped)-1].ChromEnd]...)
	}

	return ans
}

func filterZeroCoverage(fa []*fasta.Fasta, bedFile string, trim string) {
	output := fileio.EasyCreate(trim)
	defer output.Close()
	b := make(chan *bed.Bed)
	go bed.ReadToChan(bedFile, b)
	bedMap := make(map[string][]*bed.Bed)
	for i := range b {
		bedMap[i.Chrom] = append(bedMap[i.Chrom], i)
		//log.Printf("%v\n", i)
	}
	curr := &fasta.Fasta{}
	for _, chr := range fa {
		var lastPos int = 0
		curr = &fasta.Fasta{Name: chr.Name}
		for i :=0; i < len(bedMap[chr.Name]);i++ {
			//log.Printf("%s\n", dna.BasesToString(chr.Seq))
			//if dna.CountBaseInterval(chr.Seq, dna.N, int(bedMap[chr.Name][i].ChromStart), int(bedMap[chr.Name][i].ChromEnd)) == 0 {
				curr.Seq = append(curr.Seq, chr.Seq[lastPos:bedMap[chr.Name][i].ChromStart]...)
				curr.Seq = append(curr.Seq, dna.CreateAllNs(100)...)

			//}
			lastPos = int(bedMap[chr.Name][i].ChromEnd)
		}
		//append the remainder
		curr.Seq = append(curr.Seq, chr.Seq[lastPos:]...)
		//log.Printf("Seq=%s\n", dna.BasesToString(curr.Seq))
		fasta.WriteFasta(output, curr, 50)
	}



}*/
