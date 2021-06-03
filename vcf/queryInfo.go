package vcf

import (
	"github.com/vertgenlab/gonomics/exception"
	"log"
	"strconv"
	"strings"
)

// ParseInfo parses the data stored in vcf.Info.
// Fills an unexported map in the Vcf struct
// that can be queried with Query(type) functions.
func ParseInfo(v Vcf, header Header) Vcf {
	v.parsedInfo = make(map[string]interface{})
	if v.Info == "." {
		return v
	}
	var fields, tagValuePair []string
	fields = strings.Split(v.Info, ";")
	for i := range fields {
		tagValuePair = strings.Split(fields[i], "=")
		if len(tagValuePair) > 2 {
			log.Fatalf("literal '=' is illegal within an info field.\n" +
				"Found in following line:\n%s\n", v)
		}

		tag, ok := header.Info[tagValuePair[0]]
		if !ok {
			log.Fatalf("Info tag '%s' is not described in file header", tagValuePair[0])
		}

		if tag.number == "0" { // all Flag values must have number == 0
			v.parsedInfo[tagValuePair[0]] = true
			continue
		}

		v.parsedInfo[tagValuePair[0]] = parseValue(v, tagValuePair[1:], tag.Key)
	}
	return v
}

// ParseFormat parses the data stored in vcf.Format.
// Fills an unexported map in the Vcf struct
// that can be queried with Query(type) functions.
func ParseFormat(v Vcf, header Header) Vcf {
	v.parsedFormat = make(map[string]interface{})
	if len(v.Format) == 0 {
		return v
	}

	currValues := make([]string, len(v.Samples))
	for i := range v.Format {
		if v.Format[i] == "GT" { // Genotype is parsed separately
			continue
		}

		for j := range v.Samples {
			currValues[j] = v.Samples[j].FormatData[i]
		}

		tag, ok := header.Format[v.Format[i]]
		if !ok {
			log.Fatalf("Format tag '%s' is not described in file header", v.Format[i])
		}

		v.parsedFormat[v.Format[i]] = parseValue(v, currValues, tag.Key)
	}
	return v
}

// parseValue parses a []string according to the number and type defined in the header.
// each element in the slice corresponds to the value for a sample when query the Format field.
// len of s is always 1 when querying the Info field.
func parseValue(v Vcf, s []string, k Key) interface{} {
	var stringValues []string
	var err error
	switch k.dataType {
	case typeInteger:
		data := make([][]int, len(s))
		var val int
		for i := range s {
			stringValues = strings.Split(s[i], ",")
			for j := range stringValues {
				val, err = strconv.Atoi(stringValues[j])
				exception.PanicOnErr(err)
				data[i] = append(data[i], val)
			}
		}

		for i := range data {
			if !checkNumber(v, k, len(data[i])) {
				log.Panicf("unexpected number of values")
			}
		}

		return data

	case typeFloat:
		var val float64
		data := make([][]float64, len(s))
		for i := range s {
			stringValues = strings.Split(s[i], ",")
			for j := range stringValues {
				val, err = strconv.ParseFloat(stringValues[j], 64)
				exception.PanicOnErr(err)
				data[i] = append(data[i], val)
			}
		}

		for i := range data {
			if !checkNumber(v, k, len(data[i])) {
				log.Panicf("unexpected number of values")
			}
		}

		return data

	case typeString:
		data := make([][]string, len(s))
		for i := range s {
			stringValues = strings.Split(s[i], ",")
			data[i] = stringValues
		}

		for i := range data {
			if !checkNumber(v, k, len(data[i])) {
				log.Panicf("unexpected number of values")
			}
		}

		return data

	case typeCharacter:
		data := make([][]rune, len(s))
		for i := range s {
			stringValues = strings.Split(s[i], ",")
			for j := range stringValues {
				if len(stringValues[j]) > 1 {
					log.Panicf("could not convert %s to a rune\n", stringValues[j])
				}
				data[i] = append(data[i], rune(stringValues[j][0]))
			}
		}

		for i := range data {
			if !checkNumber(v, k, len(data[i])) {
				log.Panicf("unexpected number of values")
			}
		}

		return data

	default:
		log.Panicln("unknown type")
		return nil
	}
}

func checkNumber(v Vcf, k Key, length int) bool {
	switch k.number {
	case "A": // == num alt alleles
		if length != len(v.Alt) {
			return false
		}

	case "R": // == num ref + alt alleles
		if length != len(v.Alt) + 1 {
			return false
		}

	case "G": // one value for each possible genotype
	// TODO I think changes to the way GT is parsed is needed before ploidy can be determined
		return true

	case ".": // wildcard. they never make it easy do they...
		return true

	default:
		num, err := strconv.Atoi(k.number)
		if err != nil {
			log.Panicf("'%s' is not a valid number for header info", k.number)
		}
		if length != num {
			return false
		}
	}
	return true
}

// QueryInt retrieves integer values stored in the Info or Format fields of a vcf record.
// The input is a Key struct which is retrieved from the header and is keyed by Id.
// (e.g. header.Format["GT"].Key). QueryInt cannot be used until the requested field
// (Info or Format) has been parsed using the ParseInfo and ParseFormat functions.
//
// The return is a slice of slices where the first slice corresponds to the sample
// (this is always len == 1 when querying the Info field) and the second slice corresponds
// to multiple values that may be present for the given tag (e.g. ref/alt read depth may be "9,1").
func QueryInt(v Vcf, k Key) [][]int {
	if k.dataType != typeInteger {
		log.Panicf("requested QueryInt but key records data type as '%s'", k.dataType)
	}

	interfValue := query(v, k) // query value and store resulting interface
	value, ok := interfValue.([][]int) // assert value to the expected type
	if !ok { // panic if interfValue type does not match the expected type
		log.Panicf("value for tag '%s' in vcf at position '%d' " +
			"is not a [][]int as expected by header", k.Id, v.Pos)
	}
	return value
}

// QueryFloat retrieves float64 values stored in the Info or Format fields of a vcf record.
// The input is a Key struct which is retrieved from the header and is keyed by Id.
// (e.g. header.Format["GT"].Key). QueryFloat cannot be used until the requested field
// (Info or Format) has been parsed using the ParseInfo and ParseFormat functions.
//
// The return is a slice of slices where the first slice corresponds to the sample
// (this is always len == 1 when querying the Info field) and the second slice corresponds
// to multiple values that may be present for the given tag (e.g. ref/alt read depth may be "9,1").
func QueryFloat(v Vcf, k Key) [][]float64 {
	if k.dataType != typeInteger {
		log.Panicf("requested QueryFloat but key records data type as '%s'", k.dataType)
	}

	interfValue := query(v, k) // query value and store resulting interface
	value, ok := interfValue.([][]float64) // assert value to the expected type
	if !ok { // panic if interfValue type does not match the expected type
		log.Panicf("value for tag '%s' in vcf at position '%d' " +
			"is not a [][]float64 as expected by header", k.Id, v.Pos)
	}
	return value
}

// QueryFlag retrieves boolean value stored in the Info or Format fields of a vcf record.
// The input is a Key struct which is retrieved from the header and is keyed by Id.
// (e.g. header.Format["GT"].Key). QueryInt cannot be used until the requested field
// (Info or Format) has been parsed using the ParseInfo and ParseFormat functions.
//
// Note that flags are not valid in the Format field, so this query is only for Info.
func QueryFlag(v Vcf, k Key) bool {
	if k.dataType != typeInteger {
		log.Panicf("requested QueryFlag but key records data type as '%s'", k.dataType)
	}

	interfValue := query(v, k) // query value and store resulting interface
	value, ok := interfValue.(bool) // assert value to the expected type
	if !ok { // panic if interfValue type does not match the expected type
		log.Panicf("value for tag '%s' in vcf at position '%d' " +
			"is not a bool as expected by header", k.Id, v.Pos)
	}
	return value
}

// QueryString retrieves string values stored in the Info or Format fields of a vcf record.
// The input is a Key struct which is retrieved from the header and is keyed by Id.
// (e.g. header.Format["GT"].Key). QueryString cannot be used until the requested field
// (Info or Format) has been parsed using the ParseInfo and ParseFormat functions.
//
// The return is a slice of slices where the first slice corresponds to the sample
// (this is always len == 1 when querying the Info field) and the second slice corresponds
// to multiple values that may be present for the given tag (e.g. ref/alt read depth may be "9,1").
func QueryString(v Vcf, k Key) [][]string {
	if k.dataType != typeInteger {
		log.Panicf("requested QueryString but key records data type as '%s'", k.dataType)
	}

	interfValue := query(v, k) // query value and store resulting interface
	value, ok := interfValue.([][]string) // assert value to the expected type
	if !ok { // panic if interfValue type does not match the expected type
		log.Panicf("value for tag '%s' in vcf at position '%d' " +
			"is not a [][]string as expected by header", k.Id, v.Pos)
	}
	return value
}

// QueryRune retrieves rune values stored in the Info or Format fields of a vcf record.
// The input is a Key struct which is retrieved from the header and is keyed by Id.
// (e.g. header.Format["GT"].Key). QueryRune cannot be used until the requested field
// (Info or Format) has been parsed using the ParseInfo and ParseFormat functions.
//
// The return is a slice of slices where the first slice corresponds to the sample
// (this is always len == 1 when querying the Info field) and the second slice corresponds
// to multiple values that may be present for the given tag (e.g. ref/alt read depth may be "9,1").
func QueryRune(v Vcf, k Key) [][]rune {
	if k.dataType != typeInteger {
		log.Panicf("requested QueryRune but key records data type as '%s'", k.dataType)
	}

	interfValue := query(v, k) // query value and store resulting interface
	value, ok := interfValue.([][]rune) // assert value to the expected type
	if !ok { // panic if interfValue type does not match the expected type
		log.Panicf("value for tag '%s' in vcf at position '%d' " +
			"is not a [][]rune as expected by header", k.Id, v.Pos)
	}
	return value
}

// query performs the initial look and returns a value to
// a wrapper function that handles the type assertion.
func query(v Vcf, k Key) interface{} {
	var queryMap map[string]interface{}
	if k.isFormat {
		queryMap = v.parsedFormat
	} else {
		queryMap = v.parsedInfo
	}

	if queryMap == nil {
		log.Panic("requested field has not been initialized\n" +
			"The Info and Format fields must be initialized with " +
			"ParseInfo/ParseFormat before querying the respective field.\n")
	}
	return queryMap[k.Id]
}
