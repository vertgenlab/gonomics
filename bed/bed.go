package bed

import (
	"fmt"
	"os"
)

type Bed struct {
	Name  string
	Start int
	End   int
}

func WriteToFileHandle(file *os.File, records []Bed) error {
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
