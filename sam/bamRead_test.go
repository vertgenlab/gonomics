package sam

import (
	"fmt"
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

	var err error
	var curr Sam
	var actual []Sam
	for {
		curr, _, err = DecodeBam(r)
		if err == io.EOF {
			break
		}
		fmt.Sprintln(curr)
		actual = append(actual, curr)
	}

	if !equalExceptExtra(actual, expected) {
		t.Errorf("problem reading bam")
	}

	err = r.Close()
	if err != nil {
		t.Errorf(err.Error())
	}
}

func equalExceptExtra(a, b []Sam) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].QName != b[i].QName {
			return false
		}
		if a[i].Flag != b[i].Flag {
			return false
		}
		if a[i].RName != b[i].RName {
			return false
		}
		if a[i].Pos != b[i].Pos {
			return false
		}
		if a[i].MapQ != b[i].MapQ {
			return false
		}
		if cigar.ToString(a[i].Cigar) != cigar.ToString(b[i].Cigar) {
			return false
		}
		if a[i].RNext != b[i].RNext {
			return false
		}
		if a[i].PNext != b[i].PNext {
			return false
		}
		if a[i].TLen != b[i].TLen {
			return false
		}
		if dna.CompareSeqsIgnoreCase(a[i].Seq, b[i].Seq) != 0 {
			return false
		}
		if a[i].Qual != b[i].Qual {
			return false
		}
	}
	return true
}

const bigBam string = "/Users/danielsnellings/Desktop/1k.bam"

func TestBigBam(t *testing.T) {
	r, _ := OpenBam(bigBam)
	var s Sam
	var err error
	for {
		s, _, err = DecodeBam(r)
		if err == io.EOF {
			break
		}
		fmt.Fprint(ioutil.Discard, s)
	}
	r.Close()
}

func BenchmarkGonomicsBamRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r, _ := OpenBam(bigBam)
		var s Sam
		var err error
		for {
			s, _, err = DecodeBam(r)
			if err == io.EOF {
				break
			}
			fmt.Fprint(ioutil.Discard, s)
		}
		r.Close()
	}
}

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
