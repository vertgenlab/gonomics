package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
)

type Node struct {
	Id int64
	Seq  []dna.Base
}

/*
type Edge struct {
	//create later
}*/

func Read(filename string) []*Node {
	var line string
	var currSeq []dna.Base
	var answer []*Node
	var seqIdx int64 = -1
	var doneReading bool = false

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			tmp := Node{Id: common.StringToInt64(line[1:len(line)])}
			answer = append(answer, &tmp)
			seqIdx++
		} else {
			currSeq = dna.StringToBases(line)
			answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq...)
		}
	}
	return answer
}

func WriteToFileHandle(file io.Writer, records []*Node, lineLength int) {
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, ">%d\n", rec.Id)
		common.ExitIfError(err)
		for i := 0; i < len(rec.Seq); i += lineLength {
			if i+lineLength > len(rec.Seq) {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:]))
				common.ExitIfError(err)
			} else {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:i+lineLength]))
				common.ExitIfError(err)
			}
		}
	}
}

func Write(filename string, records []*Node) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteToFileHandle(file, records, lineLength)
}

