package main

import (
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/vcf"
)

func writeTableHeader(outfile io.Writer, header vcf.Header, maxAlts int) (infoOrder []vcf.InfoHeader, formatOrder []vcf.FormatHeader) {
	s := new(strings.Builder)
	if len(header.Text) == 0 {
		log.Fatal("ERROR: no vcf header found. Must have well-formed header to use -tsv")
	}

	// write common header components
	s.WriteString("Chromosome,Position,ID,Reference")
	if maxAlts == 1 {
		s.WriteString(",Alternate")
	} else {
		for i := 0; i < maxAlts; i++ {
			s.WriteString(",Alternate_" + fmt.Sprintf("%d", i))
		}
	}
	s.WriteString(",Quality,Filter")

	// write INFO field
	var numFields int
	for key, val := range header.Info {
		numFields = numberOfFields(maxAlts, val.Key)
		if numFields == 1 {
			s.WriteString("," + key)
		} else {
			for i := 0; i < numFields; i++ {
				s.WriteString("," + key + fmt.Sprintf("_%d", i))
			}
		}
		infoOrder = append(infoOrder, val)
	}

	// get formatOrder
	for _, val := range header.Format {
		formatOrder = append(formatOrder, val)
	}

	// get sample order
	sampleOrder := make([]string, len(header.Samples))
	for key, val := range header.Samples {
		sampleOrder[val] = key
	}

	// write format per sample
	for _, fmtHeader := range formatOrder {
		for _, sample := range sampleOrder {
			numFields = numberOfFields(maxAlts, fmtHeader.Key)
			if numFields == 1 {
				s.WriteString("," + fmtHeader.Id + "_" + sample)
			} else {
				for i := 0; i < numFields; i++ {
					s.WriteString("," + fmtHeader.Id + "_" + sample + fmt.Sprintf("_%d", i))
				}
			}
		}
	}

	_, err := fmt.Fprintln(outfile, s.String())
	exception.PanicOnErr(err)
	return
}

func writeAsTable(s *strings.Builder, outfile io.Writer, v vcf.Vcf, header vcf.Header, infoOrder []vcf.InfoHeader, formatOrder []vcf.FormatHeader, maxAlts int) {
	s.Reset()

	// write basic data
	s.WriteString(fmt.Sprintf("%s,%d,%s,%s,%s",
		v.Chr, v.Pos, v.Id, v.Ref, strings.Join(v.Alt, ",")))
	for i := len(v.Alt); i < maxAlts; i++ {
		s.WriteString(",")
	}
	s.WriteString(fmt.Sprintf(",%g,%s", v.Qual, v.Filter))

	// write info data
	v = vcf.ParseInfo(v, header)
	for i := range infoOrder {
		writeData(s, v, infoOrder[i].Key, numberOfFields(maxAlts, infoOrder[i].Key), 1)
	}

	// write format data
	v = vcf.ParseFormat(v, header)
	for i := range formatOrder {
		writeData(s, v, formatOrder[i].Key, numberOfFields(maxAlts, formatOrder[i].Key), len(v.Samples))
	}

	_, err := fmt.Fprintln(outfile, s.String())
	exception.PanicOnErr(err)
}

// getMaxAltCount reads through the input vcf file to determine the maximum number of alternate alleles present.
func getMaxAltCount(infile string) int {
	var maxAlts int
	records, _ := vcf.GoReadToChan(infile)
	for v := range records {
		if len(v.Alt) > maxAlts {
			maxAlts = len(v.Alt)
		}
	}
	return maxAlts
}

func numberOfFields(maxAlts int, k vcf.Key) int {
	switch k.Number {
	case "A": // == num alt alleles
		return maxAlts

	case "R": // == num ref + alt alleles
		return maxAlts + 1

	case "G": // one value for each possible genotype
		return 1 //TODO once parser is updated

	case ".": // wildcard. they never make it easy do they...
		return 1

	default:
		num, err := strconv.Atoi(k.Number)
		if err != nil {
			log.Panicf("'%s' is not a valid Number for header info", k.Number)
		}
		return num
	}
}

func writeData(s *strings.Builder, v vcf.Vcf, key vcf.Key, numberOfFieldsPerSample int, repeats int) {
	var innerFieldsWritten, fieldsWritten, i, j int

	switch key.DataType {
	case vcf.Integer:
		ints, found := vcf.QueryInt(v, key)
		if !found {
			break
		}
		for i = range ints { // per sample data
			for j = range ints[i] { // per field data
				s.WriteString(fmt.Sprintf(",%d", ints[i][j]))
				fieldsWritten++
				innerFieldsWritten++
			}
			for j = innerFieldsWritten; j < numberOfFieldsPerSample; j++ {
				s.WriteString(",")
				fieldsWritten++
			}
			innerFieldsWritten = 0
		}

	case vcf.Float:
		flts, found := vcf.QueryFloat(v, key)
		if !found {
			break
		}
		for i = range flts { // per sample data
			for j = range flts[i] { // per field data
				s.WriteString(fmt.Sprintf(",%g", flts[i][j]))
				fieldsWritten++
				innerFieldsWritten++
			}
			for j = innerFieldsWritten; j < numberOfFieldsPerSample; j++ {
				s.WriteString(",")
				fieldsWritten++
			}
			innerFieldsWritten = 0
		}

	case vcf.String:
		strs, found := vcf.QueryString(v, key)
		if !found {
			break
		}
		for i = range strs { // per sample data
			for j = range strs[i] { // per field data
				s.WriteString(fmt.Sprintf(",%s", strs[i][j]))
				fieldsWritten++
				innerFieldsWritten++
			}
			for j = innerFieldsWritten; j < numberOfFieldsPerSample; j++ {
				s.WriteString(",")
				fieldsWritten++
			}
			innerFieldsWritten = 0
		}

	case vcf.Character:
		chars, found := vcf.QueryRune(v, key)
		if !found {
			break
		}
		for i = range chars { // per sample data
			for j = range chars[i] { // per field data
				s.WriteString(fmt.Sprintf(",%c", chars[i][j]))
				fieldsWritten++
				innerFieldsWritten++
			}
			for j = innerFieldsWritten; j < numberOfFieldsPerSample; j++ {
				s.WriteString(",")
				fieldsWritten++
			}
			innerFieldsWritten = 0
		}

	case vcf.Flag:
		found := vcf.QueryFlag(v, key)
		if found {
			s.WriteString(",TRUE")
		} else {
			s.WriteString(",FALSE")
		}
		fieldsWritten++
		innerFieldsWritten++
		for j = innerFieldsWritten; j < numberOfFieldsPerSample; j++ {
			s.WriteString(",")
			fieldsWritten++
		}
		innerFieldsWritten = 0
	}

	for j = fieldsWritten; j < numberOfFieldsPerSample*repeats; j++ {
		s.WriteString(",")
	}
}
