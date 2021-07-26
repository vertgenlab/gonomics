package sam

import (
	"fmt"
	//bgBam "github.com/biogo/hts/bam"
	//bgSam "github.com/biogo/hts/sam"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"io"
	"io/ioutil"
	"os/exec"
	"testing"
)

const bamTestfile string = "../bgzf/testdata/test.bam"
const samTestfile string = "../bgzf/testdata/test.sam"

func TestReadBam(t *testing.T) {
	var err error
	r, header := OpenBam(bamTestfile)
	expected, expHeader := Read(samTestfile)

	// -1 for samtools view line that made the sam file
	if len(header.Text) != len(expHeader.Text)-1 {
		t.Errorf("problem reading bam")
	}

	for i := range header.Text {
		if header.Text[i] != expHeader.Text[i] {
			t.Errorf("problem reading bam")
		}
	}

	if len(header.Chroms) != len(r.refs) {
		t.Errorf("plain text and bam header chroms do not match")
	}

	for i := range r.refs {
		if r.refs[i] != header.Chroms[i] {
			t.Errorf("plain text and bam header chroms do not match")
		}
	}

	var actual []Sam
	for {
		var newSam Sam
		_, err = DecodeBam(r, &newSam)
		if err == io.EOF {
			break
		}
		actual = append(actual, newSam)
	}

	if !equalExceptExtra(actual, expected) {
		t.Errorf("problem reading bam")
	}

	err = r.Close()
	if err != nil {
		t.Errorf(err.Error())
	}

	if len(actual) != len(expected) {
		t.Error("problem reading bam")
	}

	for i := range actual {
		if actual[i].String() != expected[i].String() {
			t.Error("problem reading bam")
		}
	}
}

func equalExceptExtra(a, b []Sam) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].QName != b[i].QName {
			fmt.Println("qname")
			return false
		}
		if a[i].Flag != b[i].Flag {
			fmt.Println("name")
			return false
		}
		if a[i].RName != b[i].RName {
			fmt.Println("rname")
			return false
		}
		if a[i].Pos != b[i].Pos {
			fmt.Println("pos")
			return false
		}
		if a[i].MapQ != b[i].MapQ {
			fmt.Println("mapq")
			return false
		}
		if cigar.ToString(a[i].Cigar) != cigar.ToString(b[i].Cigar) {
			fmt.Println("cig")
			return false
		}
		if a[i].RNext != b[i].RNext {
			fmt.Println("rnext")
			return false
		}
		if a[i].PNext != b[i].PNext {
			fmt.Println("p")
			return false
		}
		if a[i].TLen != b[i].TLen {
			fmt.Println("t")
			return false
		}
		if dna.CompareSeqsIgnoreCase(a[i].Seq, b[i].Seq) != 0 {
			fmt.Println("seq")
			return false
		}
		if a[i].Qual != b[i].Qual {
			fmt.Println("qual")
			return false
		}
	}
	return true
}

const bigBam string = "/Users/danielsnellings/Desktop/10k.bam"
const bigSam string = "/Users/danielsnellings/Desktop/10k.sam"

func BenchmarkBamOpenClose(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r, _ := OpenBam(bigBam)
		r.Close()
	}
}

func BenchmarkBamAllocs(b *testing.B) {
	r, _ := OpenBam(bigBam)
	var s Sam
	for i := 0; i < b.N; i++ {
		_, _ = DecodeBam(r, &s)
		b.StopTimer()
		if i%9000 == 0 {
			r.Close()
			r, _ = OpenBam(bigBam)
		}
		b.StartTimer()
	}
	r.Close()
}

func BenchmarkGonomicsBamRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r, _ := OpenBam(bigBam)
		var s Sam
		var err error
		for {
			_, err = DecodeBam(r, &s)
			if err == io.EOF {
				break
			}
			fmt.Fprint(ioutil.Discard, s)
		}
		r.Close()
	}
}

//func BenchmarkBiogoBamRead(b *testing.B) {
//	for i := 0; i < b.N; i++ {
//		file, _ := os.Open(bigBam)
//		r, _ := bgBam.NewReader(file, 0)
//		var s *bgSam.Record
//		var err error
//		for {
//			s, err = r.Read()
//			if err == io.EOF {
//				break
//			}
//			fmt.Fprint(ioutil.Discard, s)
//		}
//		r.Close()
//		file.Close()
//	}
//}

func BenchmarkSamtoolsBamRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		cmd := exec.Command("samtools", "view", bigBam)
		cmd.Stdout = ioutil.Discard

		b.StartTimer()
		err := cmd.Run()
		if err != nil {
			panic(err)
		}
		b.StopTimer()
	}
}

//func BenchmarkPrint(b *testing.B) {
//	var a uint32 = 10
//	for i := 0; i < b.N; i++ {
//		print(a)
//	}
//}
//
//func BenchmarkFmtPrint(b *testing.B) {
//	var a uint32 = 10
//	for i := 0; i < b.N; i++ {
//		fmt.Print(a)
//	}
//}
