package fasta

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"os"
	"strings"
)

type Fasta struct {
	Name string
	Seq  []dna.Base
}

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
}

func ReadNew(filename string) ([]*Fasta, error) {
	var line string
	var currSeq []dna.Base
	var answer []*Fasta
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
			tmp := Fasta{Name: line[1:len(line)]}
			answer = append(answer, &tmp)
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
}
/*
func SplitRead(filename stringe, n int) ([]*Fasta, error) {
	var line string
	var currSeq []dna.Base
	var answer []*Fasta
	var seqIdx int64 = -1

	file, err := os.Open(filenam)
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
			//Append the number of Fasta slices, given we group the bases by 'n'
			for i :=0; i < len(dna.StringToBases(line))/n; i++ {
				tmp := Fasta{Name: line[1:len(line)]}
				answer = append(answer, &tmp)
				seqIdx++
			}
		default:
			currSeq, err = dna.StringToBases(line)
			if err != nil {
				log.Fatal(err)
			}
			for i:=0; i < len(currSeq)/n; i+=n {
				//Loop through and append the groups one by 1.
				answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq[i:i+n]...)
				
			}
			
		}
	}
	return answer, scanner.Err()
}
*/

func WriteToFileHandle(file *os.File, records []Fasta, lineLength int) error {
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, ">%s\n", rec.Name)
		for i := 0; i < len(rec.Seq); i += lineLength {
			if i+lineLength > len(rec.Seq) {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:]))
				if err != nil {
					return err
				}
			} else {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:i+lineLength]))
				if err != nil {
					return err
				}
			}
		}
	}
	return nil
}

func Write(filename string, records []Fasta) error {
	lineLength := 50
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	return WriteToFileHandle(file, records, lineLength)
}

func WriteGroups(filename string, groups [][]Fasta) error {
	lineLength := 50
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	for i, _ := range groups {
		err := WriteToFileHandle(file, groups[i], lineLength)
		if err != nil {
			return err
		}
		_, err = fmt.Fprint(file, "\n")
		if err != nil {
			return err
		}
	}
	return nil
}
