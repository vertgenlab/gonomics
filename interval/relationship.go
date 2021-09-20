package interval

import (
	"fmt"
	"log"
)

const (
	xMin = 0
	xMax = 10000000000 // maximum size of genome
)

func TestValidRelationship(op string) bool {
	switch op {
	case "o":
	case "oi":
	case "d":
	case "di":
	case "m":
	case "mi":
	case "s":
	case "si":
	case "f":
	case "fi":
	case "lt":
	case "gt":
	case "e":
	case "any":
	case "within":
	case "start":
	case "end":
	case "equal":
	default:
		log.Fatalf("ERROR: Invalid relationship: %s", op)
	}
	return true
}

func transform(query Interval, op string) (x1, x2, y1, y2 float64) {
	var x float64 = float64(query.GetChromStart())
	var y float64 = float64(query.GetChromEnd()-1)

	switch op {
	case "o":
		x1, y1 = xMin, x+0.5
		x2, y2 = x-0.5, y-0.5
	case "oi":
		x1, y1 = x+0.5, y+0.5
		x2, y2 = y-0.5, xMax
	case "d":
		x1, y1 = x+0.5, x+0.5
		x2, y2 = y-0.5, y-0.5
	case "di":
		x1, y1 = xMin, y+0.5
		x2, y2 = x-0.5, xMax
	case "m":
		x1, y1 = xMin, x
		x2, y2 = x, x
	case "mi":
		x1, y1 = y, y
		x2, y2 = y, xMax
	case "s":
		x1, y1 = x, x
		x2, y2 = x, y-0.5
	case "si":
		x1, y1 = x, y+0.5
		x2, y2 = x, xMax
	case "f":
		x1, y1 = x+0.5, y
		x2, y2 = y, y
	case "fi":
		x1, y1 = xMin, y
		x2, y2 = x-0.5, y
	case "lt":
		x1, y1 = xMin, xMin
		x2, y2 = x-0.5, x-0.5
	case "gt":
		x1, y1 = y+0.5, y+0.5
		x2, y2 = xMax, xMax
	case "e":
		x1, y1 = x, y
		x2, y2 = x, y
	default:
		log.Fatalf("ERROR: Invalid relationship: %s", op)
	}

	return x1, x2, y1, y2
}

func PrintRelationships() {
	o := fmt.Sprintf("\t*------*\n\t    *------*\n\n\n")
	oi := fmt.Sprintf("\t    *------*\n\t*------*\n\n\n")
	d := fmt.Sprintf("\t  *--*\n\t*------*\n\n\n")
	di := fmt.Sprintf("\t*------*\n\t  *--*\n\n\n")
	m := fmt.Sprintf("\t*------*\n\t       *------*\n\n\n")
	mi := fmt.Sprintf("\t       *------*\n\t*------*\n\n\n")
	s := fmt.Sprintf("\t*---*\n\t*------*\n\n\n")
	si := fmt.Sprintf("\t*------*\n\t*---*\n\n\n")
	f := fmt.Sprintf("\t   *---*\n\t*------*\n\n\n")
	fi := fmt.Sprintf("\t*------*\n\t   *---*\n\n\n")
	gt := fmt.Sprintf("\t*---*\n\t        *---*\n\n\n")
	lt := fmt.Sprintf("\t        *---*\n\t*---*\n\n\n")
	e := fmt.Sprintf("\t*------*\n\t*------*\n\n\n")

	fmt.Println("Valid relationships are as follows")
	fmt.Println("Top line:    target")
	fmt.Println("Bottom line: query")
	fmt.Printf("\"o\"  = %s", o)
	fmt.Printf("\"oi\" = %s", oi)
	fmt.Printf("\"d\"  = %s", d)
	fmt.Printf("\"di\" = %s", di)
	fmt.Printf("\"m\"  = %s", m)
	fmt.Printf("\"mi\" = %s", mi)
	fmt.Printf("\"s\"  = %s", s)
	fmt.Printf("\"si\" = %s", si)
	fmt.Printf("\"f\"  = %s", f)
	fmt.Printf("\"fi\" = %s", fi)
	fmt.Printf("\"gt\" = %s", gt)
	fmt.Printf("\"lt\" = %s", lt)
	fmt.Printf("\"e\"  = %s", e)

	fmt.Println("Compound relationships may be called as follows:")
	fmt.Println("\"any\"    = o + oi + d + di + m + mi + s + si + f + fi + e")
	fmt.Println("\"within\" = d + s  + f + e")
	fmt.Println("\"start\"  = s + si + e")
	fmt.Println("\"end\"    = f + fi + e")
	fmt.Println("\"equal\"  = e")
}
