package sam

import (
	"bytes"
	"encoding/hex"
	"errors"
	"fmt"
	"log"
	"math"
	"strings"
)

// QueryTag returns the value of tag t for read s. Currently QueryTag
// is only supported if the sam was read from a bam file. If no tag
// data is present, or if the record was not from a bam file, QueryTag
// will return an error. The input tag t must be exactly 2 characters.
//
// The returned value is an interface that may store any of: uint8, int8,
// uint16, int16, uint32, int32, float32, string, or a slice of interface{}
// with any of the previous underlying types.
// It is the callers responsibility to know/determine the type
// and perform the assertion.
//
// The second return is false if the requested tag was not found.
func QueryTag(s Sam, t string) (value interface{}, found bool, error error) {
	if s.parsedExtra == nil {
		s, error = parseExtra(s)
	}
	var idx int
	idx, found = s.parsedExtraIdx[Tag{t[0], t[1]}]
	if found {
		value = s.parsedExtra[idx]
	}
	return
}

// ParseExtra generates the text representation of the Extra field.
// This is required if the file was read from a bam file and the Extra
// field is going to be modified.
func ParseExtra(s *Sam) error {
	var err error
	if s.parsedExtra == nil {
		var tmp Sam
		tmp, err = parseExtra(*s)
		s.parsedExtra = tmp.parsedExtra
		s.parsedExtraIdx = tmp.parsedExtraIdx
		s.parsedExtraTags = tmp.parsedExtraTags
		s.parsedExtraTypes = tmp.parsedExtraTypes
	}
	s.Extra = parsedExtraToString(s)
	s.unparsedExtra = nil
	return err
}

func parseExtra(s Sam) (Sam, error) {
	if s.unparsedExtra == nil {
		return s, errors.New("no tags present, or record was not parsed from a bam file")
	}
	s.parsedExtraIdx = make(map[Tag]int)
	r := bytes.NewBuffer(s.unparsedExtra)

	var tag Tag
	var typ byte
	for i := 0; r.Len() > 0; i++ {
		tag[0], _ = r.ReadByte()
		tag[1], _ = r.ReadByte()
		typ, _ = r.ReadByte()
		s.parsedExtraTypes = append(s.parsedExtraTypes, []byte{typ})
		s.parsedExtraTags = append(s.parsedExtraTags, tag)
		s.parsedExtraIdx[tag] = i
		if typ == 'B' {
			typ, _ = r.ReadByte()
			s.parsedExtraTypes[i] = append(s.parsedExtraTypes[i], typ)
			count := le.Uint32(r.Next(4))
			s.parsedExtra = append(s.parsedExtra, getVals(typ, r, int(count)))
		} else {
			s.parsedExtra = append(s.parsedExtra, getVals(typ, r, 1))
		}
	}
	return s, nil
}

func getVals(typ byte, r *bytes.Buffer, count int) interface{} {
	switch typ {
	case 'A': // rune
		if count == 1 {
			return rune(r.Next(1)[0])
		} else {
			var answer []rune
			for i := 0; i < count; i++ {
				answer = append(answer, rune(r.Next(1)[0]))
			}
			return answer
		}

	case 'c': // int8
		if count == 1 {
			return int8(r.Next(1)[0])
		} else {
			answer := make([]int8, count)
			for i := 0; i < count; i++ {
				answer[i] = int8(r.Next(1)[0])
			}
			return answer
		}

	case 'C': // uint8
		if count == 1 {
			return uint8(r.Next(1)[0])
		} else {
			answer := make([]uint8, count)
			for i := 0; i < count; i++ {
				answer[i] = uint8(r.Next(1)[0])
			}
			return answer
		}

	case 's': // int16
		if count == 1 {
			return int16(le.Uint16(r.Next(2))) // le is an alias for binary.LittleEndian
		} else {
			answer := make([]int16, count)
			for i := 0; i < count; i++ {
				answer[i] = int16(le.Uint16(r.Next(2)))
			}
			return answer
		}

	case 'S': // uint16
		if count == 1 {
			return uint16(le.Uint16(r.Next(2)))
		} else {
			answer := make([]uint16, count)
			for i := 0; i < count; i++ {
				answer[i] = uint16(le.Uint16(r.Next(2)))
			}
			return answer
		}

	case 'i': // int32
		if count == 1 {
			return int32(le.Uint32(r.Next(4)))
		} else {
			answer := make([]int32, count)
			for i := 0; i < count; i++ {
				answer[i] = int32(le.Uint32(r.Next(4)))
			}
			return answer
		}

	case 'I': // uint32
		if count == 1 {
			return uint32(le.Uint32(r.Next(4)))
		} else {
			answer := make([]uint32, count)
			for i := 0; i < count; i++ {
				answer[i] = uint32(le.Uint32(r.Next(4)))
			}
			return answer
		}

	case 'f': // float32
		if count == 1 {
			return math.Float32frombits(le.Uint32(r.Next(4)))
		} else {
			answer := make([]float32, count)
			for i := 0; i < count; i++ {
				answer[i] = math.Float32frombits(le.Uint32(r.Next(4)))
			}
			return answer
		}

	case 'Z': // string
		if count == 1 {
			val, err := r.ReadString(0x0)
			if err != nil {
				log.Panic("bam may be malformed")
			}
			return trimNulOrPanic(val)
		} else {
			answer := make([]string, count)
			for i := 0; i < count; i++ {
				val, err := r.ReadString(0x0)
				if err != nil {
					log.Panic("bam may be malformed")
				}
				answer[i] = trimNulOrPanic(val)
			}
			return answer
		}

	case 'H': // []byte
		if count == 1 {
			val, err := r.ReadBytes(0x0)
			if err != nil {
				log.Panic("bam may be malformed")
			}
			return val[:len(val)-1] // trim NUL byte
		} else {
			answer := make([][]byte, count)
			for i := 0; i < count; i++ {
				val, err := r.ReadString(0x0)
				if err != nil {
					log.Panic("bam may be malformed")
				}
				answer[i], err = hex.DecodeString(trimNulOrPanic(val))
				if err != nil {
					log.Panic("bam may be malformed")
				}
			}
			return answer
		}

	default:
		log.Panicf("unrecognized value type in bam file '%s'", string(typ))
		return nil
	}
}

func parsedExtraToString(r *Sam) string {
	var s strings.Builder
	var vals []string
	for i := range r.parsedExtraTypes {
		if r.parsedExtraTypes[i][0] == 'B' {
			switch r.parsedExtraTypes[i][1] {
			case 'c':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]int8) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'C':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]uint8) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 's':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]int16) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'S':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]uint16) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'i':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]int32) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'I':
				s.WriteString(fmt.Sprintf("%s:B:i:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]uint32) {
					vals = append(vals, fmt.Sprintf("%d", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'f':
				s.WriteString(fmt.Sprintf("%s:B:f:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([]float32) {
					vals = append(vals, fmt.Sprintf("%f", val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")

			case 'Z':
				s.WriteString(fmt.Sprintf("%s:B:Z:%s\t", r.parsedExtraTags[i], strings.Join(r.parsedExtra[i].([]string), ",")))

			case 'H':
				s.WriteString(fmt.Sprintf("%s:B:H:", r.parsedExtraTags[i]))
				vals = vals[:0]
				for _, val := range r.parsedExtra[i].([][]byte) {
					vals = append(vals, hex.EncodeToString(val))
				}
				s.WriteString(strings.Join(vals, ",") + "\t")
			}
		} else {
			switch r.parsedExtraTypes[i][0] {
			case 'c':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(int8)))

			case 'C':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(uint8)))

			case 's':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(int16)))

			case 'S':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(uint16)))

			case 'i':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(int32)))

			case 'I':
				s.WriteString(fmt.Sprintf("%s:i:%d\t", r.parsedExtraTags[i], r.parsedExtra[i].(uint32)))

			case 'f':
				s.WriteString(fmt.Sprintf("%s:f:%f\t", r.parsedExtraTags[i], r.parsedExtra[i].(float32)))

			case 'Z':
				s.WriteString(fmt.Sprintf("%s:Z:%s\t", r.parsedExtraTags[i], r.parsedExtra[i].(string)))

			case 'H':
				s.WriteString(fmt.Sprintf("%s:H:%s\t", r.parsedExtraTags[i], hex.EncodeToString(r.parsedExtra[i].([]byte))))
			}
		}
	}
	return strings.Trim(s.String(), "\t")
}
