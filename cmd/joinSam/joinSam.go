package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func uni(a sam.Sam) bool {
	_, found, _ := sam.QueryTag(a, "XA")
	if !sam.IsUnmapped(a) && a.MapQ > 24 && !found {
		return true
	} else {
		return false
	}
}

func lowQualMulti(a sam.Sam) bool {
	_, found, _ := sam.QueryTag(a, "XA")
	if !sam.IsUnmapped(a) && a.MapQ < 25 && !found {
		return true
	} else {
		return false
	}
}

func joinSam(in1, in2, out, sumFile string) {
	var j sam.Sam
	var err error
	var total, unmapped, mapped, uniCount, multiLowQ, multi, singleton int
	chan1, head1 := sam.GoReadToChan(in1)
	chan2, _ := sam.GoReadToChan(in2)

	o := fileio.EasyCreate(out)
	bw := sam.NewBamWriter(o, head1)

	sumOut := fileio.EasyCreate(sumFile)

	for i := range chan1 {
		total++
		j = <-chan2
		if i.QName != j.QName {
			log.Fatalf("Alignment files are either not sorted by read-name or malformed. Read names: %s\t%s\n", i.QName, j.QName)
		}
		if sam.IsUnmapped(i) && sam.IsUnmapped(j) {
			unmapped++
			continue
		} else if !sam.IsUnmapped(i) && !sam.IsUnmapped(j) {
			mapped++
			if uni(i) && uni(j) {
				uniCount++
				i, j = updateSam(i, j)
				sam.WriteToBamFileHandle(bw, i, 0)
				sam.WriteToBamFileHandle(bw, j, 0)
				continue
			}
			if lowQualMulti(i) || lowQualMulti(j) {
				multiLowQ++
				continue
			}
			i, j = updateSam(i, j)
			sam.WriteToBamFileHandle(bw, i, 0)
			sam.WriteToBamFileHandle(bw, j, 0)
			multi++
		} else {
			singleton++
		}
	}
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Total read pairs: %d", total))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Unmapped read pairs: %d", unmapped))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Singleton read pairs: %d", singleton))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Total number of mapped reads: %d", mapped))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Uniquely mapping read pairs: %d", uniCount))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Total mulit-mapping read pairs: %d", multi+multiLowQ))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Multi-mapping, high-quality read pairs: %d", multi))
	fileio.WriteToFileHandle(sumOut, fmt.Sprintf("Multi-mapping, low-quality read pairs: %d", multiLowQ))

	err = bw.Close()
	exception.PanicOnErr(err)
	err = o.Close()
	exception.PanicOnErr(err)
	err = sumOut.Close()
	exception.PanicOnErr(err)
}

// updateSam takes in two paired sam entries and updates them to paired end
func updateSam(r1, r2 sam.Sam) (sam.Sam, sam.Sam) {
	//store flags as int
	f1 := r1.Flag
	f2 := r2.Flag
	//check and update flags
	switch {
	case f1&0x4 == 0x4: //4 = unmapped
		f1 = f1 | 0x8 //8 = unmapped but paired
	case f2&0x4 == 0x4:
		f2 = f2 | 0x8
	default: // update both reads to show that they are proper pairs
		f1 = f1 | 0x1
		f1 = f1 | 0x2
		f2 = f2 | 0x1
		f2 = f2 | 0x2
	}
	//update strand info
	if f1&0x10 == 0x10 {
		f2 = f2 | 0x20
	} else if f2&0x10 == 0x10 {
		f1 = f1 | 0x20
	}

	//first vs second read in pair
	f1 = f1 | 0x40
	f2 = f2 | 0x80

	//insert updated flags back into sam entry
	r1.Flag = f1
	r2.Flag = f2

	//add RNEXT / PNEXT info
	r1.RNext = r2.RName
	r2.RNext = r1.RName
	r1.PNext = r2.Pos
	r2.PNext = r1.Pos

	return r1, r2
}

func main() {
	flag.Parse()

	joinSam(flag.Arg(0), flag.Arg(1), flag.Arg(2), flag.Arg(3))
}
