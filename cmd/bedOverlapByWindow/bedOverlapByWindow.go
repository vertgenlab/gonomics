package main

import (
        "github.com/vertgenlab/gonomics/chromInfo"
        "github.com/vertgenlab/gonomics/bed"
        "github.com/vertgenlab/gonomics/fileio"
        "github.com/vertgenlab/gonomics/common"
        "fmt"
        "flag"
        "log"
)

func bedOverlapByWindow (infile string, chromsizes string, outfile string, windowSize int64) {
        cInfo := chromInfo.ReadToSlice(chromsizes)
        bInfo := bed.Read(infile)
        out := fileio.EasyCreate(outfile)
        defer out.Close()
        var positionCounts map[string][]uint32
        positionCounts = make(map[string][]uint32)
        var i, b, p, j, x int64
        var thisChrom []uint32 


        for i = 0; i < int64(len(cInfo)); i++ {

                positionCounts[cInfo[i].Name] = make([]uint32, cInfo[i].Size)


        }

        for b = 0; b < int64(len(bInfo)); b++ {
                
                thisChrom = positionCounts[bInfo[b].Chrom]


                for p = bInfo[b].ChromStart; p < bInfo[b].ChromEnd; p++{

                        for x = common.MaxInt64(0, p-(windowSize-1)); x < common.MinInt64(bInfo[b].ChromEnd, p+1); x++{

//for x = common.MaxInt64(0, p-(windowSize-1)); x < common.MinInt64(bInfo[b].ChromEnd, p+(windowSize-1)); x++{
                           thisChrom[x]+=1  

                        }

                }

        }
        
        for i = 0; i < int64(len(cInfo)); i++ {

                thisChrom=positionCounts[cInfo[i].Name]

                for j = 0; j < int64(len(thisChrom)); j++{

                        fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%d\n", cInfo[i].Name, j, j+windowSize, ".", thisChrom[j])

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

