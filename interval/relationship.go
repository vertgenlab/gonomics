package interval

import "log"

const (
	xMin = 0
	xMax = 10000000000 // maximum size of genome
)

func transform(query Interval, op string) (x1, x2, y1, y2 float32) {
	var x float32 = float32(query.GetChromStart())
	var y float32 = float32(query.GetChromEnd())

	switch op {
	case "o":
		x1, x2 = xMin, x + 0.5
		y1, y2 = x - 0.5, y - 0.5
	case "oi":
		x1, x2 = x + 0.5, y + 0.5
		y1, y2 = y - 0.5, xMax
	case "d":
		x1, x2 = x + 0.5, x + 0.5
		y1, y2 = y - 0.5, y - 0.5
	case "di":
		x1, x2 = xMin, y + 0.5
		y1, y2 = x - 0.5, xMax
	case "m":
		x1, x2 = xMin, x
		y1, y2 = x, x
	case "mi":
		x1, x2 = y, y
		y1, y2 = y, xMax
	case "s":
		x1, x2 = x, x
		y1, y2 = x, y - 0.5
	case "si":
		x1, x2 = x, y + 0.5
		y1, y2 = x, xMax
	case "f":
		x1, x2 = x + 0.5, y
		y1, y2 = y, y
	case "fi":
		x1, x2 = xMin, y
		y1, y2 = x - 0.5, y
	case "lt":
		x1, x2 = xMin, xMin
		y1, y2 = x - 0.5, x - 0.5
	case "gt":
		x1, x2 = y + 0.5, y + 0.5
		y1, y2 = xMax, xMax
	case "e":
		x1, x2 = x, y
		y1, y2 = x, y
	default:
		log.Fatalf("ERROR: Invalid relationship: %s", op)
	}
	return x1, x2, y1, y2
}
