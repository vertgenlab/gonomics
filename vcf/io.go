package vcf

import (
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/exception"

	"github.com/vertgenlab/gonomics/fileio"
)

// Read parses a slice of VCF structs from an input filename. Does not store or return the header.
func Read(filename string) ([]Vcf, Header) {
	var answer []Vcf
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	var curr Vcf
	var done bool
	for curr, done = NextVcf(file); !done; curr, done = NextVcf(file) {
		answer = append(answer, curr)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	return answer, header
}

// ReadToChan is a helper function of GoReadToChan.
func ReadToChan(file *fileio.EasyReader, data chan<- Vcf, wg *sync.WaitGroup) {
	for curr, done := NextVcf(file); !done; curr, done = NextVcf(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan parses VCF structs from an input filename and returns a chan of VCF structs along with the header of the VCF file.
func GoReadToChan(filename string) (<-chan Vcf, Header) {
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	var wg sync.WaitGroup
	data := make(chan Vcf, 1000)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, header
}

// processVcfLine is a helper function of NextVcf that parses a VCF struct from an input line of a VCF file provided as a string.
func processVcfLine(line string) Vcf {
	var curr Vcf
	var err error
	data := strings.Split(line, "\t")
	if len(data) < 8 {
		log.Panicf("Error: when reading this vcf line:\n%s\nExpecting at least 8 columns", line)
	}

	curr.Chr = data[0]
	curr.Pos, err = strconv.Atoi(data[1])
	if err != nil {
		log.Panicf("Error: VCF reading\nCould not convert '%s' to an integer in the following line\n%s\n", data[1], line)
	}
	curr.Id = data[2]
	curr.Ref = data[3]
	curr.Alt = strings.Split(data[4], ",")
	curr.Qual = 255
	if data[5] != "." {
		curr.Qual, err = strconv.ParseFloat(data[5], 64)
		if err != nil {
			log.Panicf("Error: VCF reading\nCould not convert '%s' to a float in the following line\n%s\n", data[5], line)
		}
	}
	curr.Filter = data[6]
	curr.Info = data[7]

	if len(data) < 9 {
		// return if no format field found
		return curr
	}
	curr.Format = strings.Split(data[8], ":")
	curr.Samples = parseSamples(data[9:], curr.Format, line)
	return curr
}

// parseSamples is a helper function of processVcfLine. Generates a slice of Sample structs from a VCF data line.
func parseSamples(samples []string, format []string, line string) []Sample {
	//DEBUG: fmt.Printf("Format: %s. Format[0]: %s.\n", format, format[0])
	if format[0] == "." || len(format) == 0 {
		return nil
	}

	for i := 1; i < len(format); i++ {
		if format[i] == "GT" {
			log.Fatalf("ERROR: The format key GT is reserved for genotype information. VCF specifications require GT to be present as"+
				"the first format key if it is present. GT key is misplaced in the following line:\n%s\n", line)
		}
	}

	answer := make([]Sample, len(samples))
	for i := range samples {
		answer[i].FormatData = strings.Split(samples[i], ":")
		if format[0] == "GT" {
			answer[i].Alleles, answer[i].Phase = parseGenotype(answer[i].FormatData[0], line)
			answer[i].FormatData[0] = ""
		}
	}
	return answer
}

// parseGenotype returns the alleles and phase parsed from the GT field in Samples.
func parseGenotype(gt string, line string) (alleles []int16, phase []bool) {
	var alleleId int64
	var err error
	// TODO: Consider reducing allocation of string.Builder variable to further improve performance.
	if gt == "." || gt == "./." {
		return nil, nil
	}

	// split GT by '/' or '|'
	text := splitGenotype(gt)
	if len(text) == 0 {
		return
	}

	alleles = make([]int16, 0, (len(text)+1)/2)
	phase = make([]bool, 1, cap(alleles)) // phase starts at size 1 for first alleles phase which is set at the end

	for i := range text {
		if i%2 == 0 { // is allele id
			if text[i] == "." {
				alleles = append(alleles, -1)
				continue
			}
			alleleId, err = strconv.ParseInt(text[i], 10, 16)
			if err != nil {
				log.Fatalf("Error: VCF reading\nCould not convert '%s' to an int16 in the following line\n%s\n", text[i], line)
			}
			alleles = append(alleles, int16(alleleId))
		} else { // is phase info
			phase = append(phase, text[i] == "|")
		}
	}

	// check to see if all the alleles so far are phased.
	// if all are phased then alleles[0] is also phased.
	// if any are false, then alleles[0] is not phased.
	var allPhased bool = true
	for i := range phase {
		if !phase[i] {
			allPhased = false
			break
		}
	}
	phase[0] = allPhased
	return
}

// splitGenotype splits each elements of the GT field into a slice of elements (e.g. 1/1 becomes []string{"1", "/", "1").
func splitGenotype(gt string) []string {
	var answer []string
	var builder strings.Builder

	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			answer = append(answer, builder.String())
			answer = append(answer, string(gt[i]))
			builder.Reset()
		} else {
			builder.WriteByte(gt[i])
		}
	}
	answer = append(answer, builder.String())
	return answer
}

// NextVcf is a helper function of Read and GoReadToChan. Checks a reader for additional data lines and parses a Vcf line if more lines exist.
func NextVcf(reader *fileio.EasyReader) (Vcf, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return Vcf{}, true
	}
	return processVcfLine(line), false
}

// TODO(craiglowe): Look into unifying WriteVcfToFileHandle and WriteVcf and benchmark speed. geno bool variable determines whether to print notes or genotypes.
func WriteVcfToFileHandle(file io.Writer, input []Vcf) {
	for i := 0; i < len(input); i++ {
		WriteVcf(file, input[i])
	}
}

// WriteVcf writes an individual Vcf struct to an io.Writer.
func WriteVcf(writer io.Writer, record Vcf) {
	_, err := fmt.Fprintf(writer, "%s\n", record.String())
	exception.PanicOnErr(err)
}

// Write writes a []Vcf to an output filename.
func Write(filename string, data []Vcf) {
	file := fileio.EasyCreate(filename)
	WriteVcfToFileHandle(file, data)
	exception.PanicOnErr(file.Close())
}

// IsVcfFile checks suffix of filename to confirm if the file is a vcf formatted file.
func IsVcfFile(filename string) bool {
	if strings.HasSuffix(filename, ".vcf") || strings.HasSuffix(filename, ".vcf.gz") {
		return true
	} else {
		return false
	}
}
