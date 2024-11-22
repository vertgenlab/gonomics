package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"os"
	"strings"
	"testing"
)

func TestTags(t *testing.T) {
	var err error
	data, _ := Read("testdata/unmapped.bam")
	for i := range data {
		err = ParseExtra(&data[i])
		if err != nil {
			t.Error("problem parsing tags in ParseExtra")
		}
	}

	expected, _ := Read("testdata/unmapped.sam")

	for i := range data {
		if data[i].Extra != expected[i].Extra {
			fmt.Println(data[i].Extra)
			fmt.Println(expected[i].Extra)
			t.Error("problem parsing tags in ParseExtra")
		}
	}
}

func TestBamWriteTags(t *testing.T) {
	file, err := os.CreateTemp("", "")
	data, header := Read("testdata/unmapped.sam")
	w := NewBamWriter(file, header)
	for i := range data {
		WriteToBamFileHandle(w, data[i], 0)
	}
	err = w.Close()
	exception.PanicOnErr(err)
	file.Close()

	var afterWrite []Sam
	br, _ := OpenBam(file.Name())
	for err != io.EOF {
		var curr Sam
		_, err = DecodeBam(br, &curr)
		ParseExtra(&curr)
		if err == nil {
			afterWrite = append(afterWrite, curr)
		}
	}
	err = br.Close()
	exception.PanicOnErr(err)

	for i := range data {
		if data[i].Extra != afterWrite[i].Extra {
			fmt.Println(data[i].Extra)
			fmt.Println(afterWrite[i].Extra)
			t.Error("problem writing tags in bam files")
		}
	}
}

func TestRemoveAllTags(t *testing.T) {
	file, err := os.CreateTemp("", "")
	data, header := Read("testdata/peak.bam")
	w := NewBamWriter(file, header)
	for i := range data {
		RemoveAllTags(&data[i])
		WriteToBamFileHandle(w, data[i], 0)
	}
	err = w.Close()
	exception.PanicOnErr(err)
	file.Close()

	var afterWrite []Sam
	br, _ := OpenBam(file.Name())
	for err != io.EOF {
		var curr Sam
		_, err = DecodeBam(br, &curr)
		ParseExtra(&curr)
		if err == nil {
			afterWrite = append(afterWrite, curr)
		}
	}
	err = br.Close()
	exception.PanicOnErr(err)

	for i := range afterWrite {
		if afterWrite[i].Extra != "" {
			t.Error("problem removing tags")
		}
	}
}

func TestRemoveTag(t *testing.T) {
	file, err := os.CreateTemp("", "")
	data, header := Read("testdata/peak.bam")
	w := NewBamWriter(file, header)
	for i := range data {
		exception.PanicOnErr(RemoveTag(&data[i], "MD")) // first tag
		exception.PanicOnErr(RemoveTag(&data[i], "RG")) // middle tag
		exception.PanicOnErr(RemoveTag(&data[i], "YT")) // last tag
		exception.PanicOnErr(RemoveTag(&data[i], "pq")) // missing tag
		WriteToBamFileHandle(w, data[i], 0)
	}
	err = w.Close()
	exception.PanicOnErr(err)
	err = file.Close()
	exception.PanicOnErr(err)

	var afterWrite []Sam
	br, _ := OpenBam(file.Name())
	for err != io.EOF {
		var curr Sam
		_, err = DecodeBam(br, &curr)
		ParseExtra(&curr)
		if err == nil {
			afterWrite = append(afterWrite, curr)
		}
	}
	err = br.Close()
	exception.PanicOnErr(err)

	if len(afterWrite) != len(data) {
		t.Error("bam files not the same length, problem writing")
	}
	for i := range afterWrite {
		if afterWrite[i].Extra != slowSafeRemoveTags(&data[i], []string{"MD", "RG", "YT"}) {
			t.Errorf("problem removing specific bam tags\nreceived:\t%s\ncorrect:\t%s", afterWrite[i].Extra, slowSafeRemoveTags(&data[i], []string{"MD", "RG", "YT"}))
			return
		}
	}
}

func TestAddTag(t *testing.T) {
	file, err := os.CreateTemp("", "")
	data, header := Read("testdata/peak.bam")
	w := NewBamWriter(file, header)
	for i := range data {
		exception.PanicOnErr(AddTag(&data[i], "pq", "Z", "<---look at those arms!"))
		WriteToBamFileHandle(w, data[i], 0)
	}
	err = w.Close()
	exception.PanicOnErr(err)
	err = file.Close()
	exception.PanicOnErr(err)

	var afterWrite []Sam
	br, _ := OpenBam(file.Name())
	for err != io.EOF {
		var curr Sam
		_, err = DecodeBam(br, &curr)
		ParseExtra(&curr)
		if err == nil {
			afterWrite = append(afterWrite, curr)
		}
	}
	err = br.Close()
	exception.PanicOnErr(err)

	if len(afterWrite) != len(data) {
		t.Error("bam files not the same length, problem writing")
	}
	for i := range afterWrite {
		if afterWrite[i].Extra != data[i].Extra {
			t.Errorf("problem removing specific bam tags\nreceived:\t%s\ncorrect:\t%s", afterWrite[i].Extra, data[i].Extra)
			return
		}
	}
}

func slowSafeRemoveTags(rec *Sam, tags []string) string {
	var err error

	// parse the extra field where tags are since we are coming from bam
	err = ParseExtra(rec)
	exception.PanicOnErr(err)

	var ans string = rec.Extra

	for j := range tags {
		splitExtra := strings.Split(ans, "\t")
		var newExtra []string = make([]string, 0, len(splitExtra)+1)

		// first remove any exiting read group and pull the tag value at the same time
		for i := range splitExtra {
			if !strings.HasPrefix(splitExtra[i], tags[j]) {
				newExtra = append(newExtra, splitExtra[i])
			}
		}
		ans = strings.Join(newExtra, "\t")
	}

	return ans
}
