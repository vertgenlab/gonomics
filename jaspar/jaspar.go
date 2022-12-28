package jaspar

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

//Pfm is a struct encoding a position frequency matrix.
type Pfm struct {
	Id   string
	Name string
	Mat  [][]float64
}

//Ppm is a struct encoding a position probability matrix.
// While we could use one struct for both objects, I think two different structs
// will help us avoid silent errors downstream.
type Ppm struct {
	Id   string
	Name string
	Mat  [][]float64
}

// WritePfmSlice writes a slice of Pfm structs to an output filename.
func WritePfmSlice(filename string, records []Pfm) {
	var err error

	file := fileio.EasyCreate(filename)

	for _, pfm := range records {
		WritePfm(file, pfm)
	}

	err = file.Close()
	exception.PanicOnErr(err)
}

// WritePfm writes an individual Pfm struct to a fileio.EasyWriter.
func WritePfm(file *fileio.EasyWriter, pfm Pfm) {
	var err error
	_, err = fmt.Fprintf(file, ">%s\t%s\n%s", pfm.Id, pfm.Name, matToString(pfm.Mat))
	exception.PanicOnErr(err)
}

// WritePpmSlice writes a slice of Ppm structs to an output filename.
func WritePpmSlice(filename string, records []Ppm) {
	var err error
	file := fileio.EasyCreate(filename)
	for _, ppm := range records {
		WritePpm(file, ppm)
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// WritePpm writes an individual Ppm struct to a fileio.EasyWriter.
func WritePpm(file *fileio.EasyWriter, ppm Ppm) {
	var err error
	_, err = fmt.Fprintf(file, ">%s\t%s\n%s", ppm.Id, ppm.Name, matToString(ppm.Mat))
	exception.PanicOnErr(err)
}

//matToString is a helper function of matToString that converts a Pfm/Ppm matrix to a string.
func matToString(mat [][]float64) string {
	var answer string = "A\t[\t"
	if len(mat) != 4 {
		log.Fatalf("Error: Input PFM must have 4 rows, one for each nucleotide.")
	}
	for i := range mat[0] {
		answer = answer + fmt.Sprintf("\t%g", mat[0][i])
	}
	answer = answer + "\t]\nC [ "
	for i := range mat[1] {
		answer = answer + fmt.Sprintf("\t%g", mat[1][i])
	}
	answer = answer + "\t]\nG [ "
	for i := range mat[2] {
		answer = answer + fmt.Sprintf("\t%g", mat[2][i])
	}
	answer = answer + "\t]\nT [ "
	for i := range mat[3] {
		answer = answer + fmt.Sprintf("\t%g", mat[3][i])
	}
	answer = answer + "\t]\n"
	return answer
}

//ReadPfm parses a slice of Pfm structs from an input file in JASPAR PFM format.
func ReadPfm(filename string) []Pfm {
	var err error
	var curr Pfm
	var answer []Pfm
	var doneReading bool
	usedMotifIds := make(map[string]bool)

	file := fileio.EasyOpen(filename)

	for curr, doneReading = NextPfm(file); !doneReading; curr, doneReading = NextPfm(file) {
		if usedMotifIds[curr.Id] {
			log.Fatalf("Error: %s is used as the ID for multiple records. IDs must be unique.", curr.Id)
		} else {
			usedMotifIds[curr.Id] = true
			answer = append(answer, curr)
		}
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// NextPfm reads and parses a single Pfm record from an input EasyReader. Returns true when the file is fully read.
func NextPfm(file *fileio.EasyReader) (Pfm, bool) {
	var header string
	var fields []string
	var motifLen int
	var answer Pfm

	line1, done1 := fileio.EasyNextRealLine(file)
	line2, done2 := fileio.EasyNextRealLine(file)
	line3, done3 := fileio.EasyNextRealLine(file)
	line4, done4 := fileio.EasyNextRealLine(file)
	line5, done5 := fileio.EasyNextRealLine(file)

	if done1 {
		return Pfm{}, true
	}
	if done2 || done3 || done4 || done5 {
		log.Fatalf("Error: There is an empty line in this Pfm record.")
	}

	if !strings.HasPrefix(line1, ">") {
		log.Fatalf("Error: Pfm header line must begin with a '>' symbol.")
	}

	//parse header line
	header = line1[1:]              //trim '>' symbol
	fields = strings.Fields(header) //split fields by whitespace ' ' or '\t', etc.
	if len(fields) == 0 {
		log.Fatalf("Error: Pfm has empty header.")
	}
	answer = Pfm{Id: fields[0]}
	if len(fields) > 1 {
		answer.Name = fields[1]
	}

	answer.Mat = make([][]float64, 4) //make 4 matrix rows, one for each nucleotide.
	motifLen = getMotifLen(line2)
	if motifLen < 1 {
		log.Fatalf("Error: Motif length must be at least 1.")
	}

	parseMotifLine(answer, line2, motifLen, 0)
	parseMotifLine(answer, line3, motifLen, 1)
	parseMotifLine(answer, line4, motifLen, 2)
	parseMotifLine(answer, line5, motifLen, 3)
	return answer, false
}

// getMotifLen takes in a line from a JASPAR PFM file and determines the length of the motif.
func getMotifLen(line string) int {
	line = strings.Replace(line, "[", " ", 1)
	line = strings.Replace(line, "]", " ", 1)
	fields := strings.Fields(line)
	return len(fields) - 1
}

// parseMotifLine is a helper function of ReadPfm that parses a line of a PFM file into.
// a line of the Pfm.Mat data structure.
func parseMotifLine(answer Pfm, line string, motifLen int, index int) {
	line = strings.Replace(line, "[", " ", 1)
	line = strings.Replace(line, "]", "", 1)
	fields := strings.Fields(line)
	if len(fields)-1 != motifLen {
		fmt.Printf("LineMotifLength: %v. MotifLen: %v.\n", len(fields)-1, motifLen)
		log.Fatalf("Error: motif length is not equal to other lines in pfm record.")
	}
	answer.Mat[index] = make([]float64, motifLen)
	fields = fields[1:] //trim first field, which corresponds to nucleotide id.
	for i := 0; i < len(fields); i++ {
		answer.Mat[index][i] = common.StringToFloat64(fields[i])
	}
}
