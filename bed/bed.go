package bed

import ("fmt"
	"os"
)

type Bed struct {
	Name string
	Start int
	End int
	
}
/*
func Read(filename string) ([]Fasta, error) {
	var line string
	var currSeq []dna.Base
	var answer []Fasta
	var seqIdx int64 = -1

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		switch {
		case strings.HasPrefix(line, "#"):
			// comment line in fasta file
		case strings.HasPrefix(line, ">"):
			answer = append(answer, Fasta{Name: line[1:len(line)]})
			seqIdx++
		default:
			currSeq, err = dna.StringToBases(line)
			if err != nil {
				log.Fatal(err)
			}
			answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq...)
		}
	}
	return answer, scanner.Err()
}*/

func WriteToFileHandle(file *os.File, records []Bed) error{
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, ">%s\t%d\t%d\n", rec.Name, rec.Start, rec.End)
		if err != nil {
			return err
		}
	}
	return nil
}

func Write(filename string, records []Bed) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	return WriteToFileHandle(file, records)
}

