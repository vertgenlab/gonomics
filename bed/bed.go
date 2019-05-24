package bed

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type Bed struct {
	Chrom      string
	ChromStart int64
	ChromEnd   int64
	Name       string
	Score      int64
	Strand     rune
}

func BedToString(bunk *Bed, fields int) string {
	switch fields {
	case 3:
		return fmt.Sprintf("%s\t%d\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd)
	case 4:
		return fmt.Sprintf("%s\t%d\t%d\t%s", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name)
	case 5:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score)
	case 6:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score, bunk.Strand)
	default:
		log.Fatalf("Error: expecting a request to print 3 to 6 bed fields, but got: %d\n", fields)
	}
	return ""
}

func WriteToFileHandle(file *os.File, records []*Bed, fields int) error {
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, "%s\n", BedToString(rec, fields))
		if err != nil {
			return err
		}
	}
	return nil
}

func Write(filename string, records []*Bed, fields int) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	return WriteToFileHandle(file, records, fields)
}

func Read(filename string) ([]*Bed, error) {
	var line string
	var answer []*Bed
	var err2 error
	var startNum, endNum int64

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line = scanner.Text()
		words := strings.Split(line, "\t")

		if len(words) < 3 || len(words) > 6 {
			return nil, fmt.Errorf("All lines in bed file should have 3 to 6 columns, but only %d were found in %s", len(words), line)
		}
		startNum, err2 = strconv.ParseInt(words[1], 10, 64)
		if err != nil {
			return nil, err
		}
		endNum, err2 = strconv.ParseInt(words[2], 10, 64)
		if err2 != nil {
			return nil, err2
		}

		current := Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
		if len(words) >= 4 {
			current.Name = words[3]
		}
		if len(words) >= 5 {
			current.Score, err2 = strconv.ParseInt(words[4], 10, 64)
			if err2 != nil {
				return nil, err2
			}
		}
		if len(words) >= 6 {
			switch words[5] {
			case "+":
				current.Strand = '+'
			case "-":
				current.Strand = '-'
			default:
				return nil, fmt.Errorf("Expecting + or - for the strand, but got %s as part of line %s", words[5], line)
			}
		}
		answer = append(answer, &current)

	}
	return answer, scanner.Err()
}
