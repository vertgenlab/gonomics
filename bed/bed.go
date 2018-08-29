package bed

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Bed struct {
	Name  string
	Start int
	End   int
}

func WriteToFileHandle(file *os.File, records []*Bed) error {
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, "%s\t%d\t%d\n", rec.Name, rec.Start, rec.End)
		if err != nil {
			return err
		}
	}
	return nil
}

func Write(filename string, records []*Bed) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	return WriteToFileHandle(file, records)
}

func Read(filename string) ([]*Bed, error) {
	var line string
	var answer []*Bed
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		words := strings.Split(line, "\t")
		startNum, err := strconv.Atoi(words[1])
		if err != nil {
			return nil, err
		}
		endNum, err2 := strconv.Atoi(words[2])
		if err2 != nil {
			return nil, err2
		}
		tmp := Bed{Name: words[0], Start: startNum, End: endNum}

		answer = append(answer, &tmp)

	}
	return answer, scanner.Err()
}
