package main

import (
        "github.com/vertgenlab/gonomics/chromInfo"
        "github.com/vertgenlab/gonomics/bed"
        "github.com/vertgenlab/gonomics/fileio"
        "fmt"
        "flag"
        "log"
)

func bedOverlapByWindow (infile string, chromsizes string, outfile string, windowSize int64) {
        cInfo := chromInfo.ReadToSlice(chromsizes)
        bInfo := bed.Read(infile)
        out := fileio.EasyCreate(outfile)
        defer out.Close()

        for i := 0; i < len(cInfo); i++ {
                for s := int64(0); s < cInfo[i].Size-windowSize; s++ {

                        output := bed.Bed{Chrom: cInfo[i].Name, ChromStart:int64(s), ChromEnd:int64(s+windowSize), Score: int64(0), Name: "."}


                        for b := 0; b < len(bInfo); b++ {
                                output.Score = output.Score+bed.OverlapLength(bInfo[b], &output)

                        }
                        bed.WriteBed(out.File, &output, 5)
                }

        }

}


func usage() {
        fmt.Print(
                "bedOverlapByWindow takes a sorted bed and counts bp in bed regions within a window size. Default is 5000bp\n" +
                        "Usage:\n" +
                        "bedOverlapByWindow input.bed chrom.sizes output.bed\n" +
                        "options:\n")
        flag.PrintDefaults()
}

func main() {
        var expectedNumArgs = 3
        var windowSize *int64 = flag.Int64("windowSize", 5000, "Specifies window size for counting bp in bed regions.")

        flag.Usage = usage
        log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
        flag.Parse()

        if len(flag.Args()) != expectedNumArgs {
                flag.Usage()
                log.Fatalf("Error: expecting %d arguments, but got %d\n",
                        expectedNumArgs, len(flag.Args()))
        }

        infile := flag.Arg(0)
        chromsizes := flag.Arg(1)
        outfile := flag.Arg(2)

        bedOverlapByWindow(infile, chromsizes, outfile, *windowSize)
}

