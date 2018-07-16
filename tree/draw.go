package tree

import (
	"fmt"
	"github.com/craiglowe/gonomics/sketch"
	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"
	"image"
	"image/color"
	//"golang.org/x/image/font/inconsolata"
)

func addLabel(img *image.RGBA, label string, xStart int, yStart int) {
	point := fixed.Point26_6{fixed.Int26_6(xStart * 64), fixed.Int26_6(yStart * 64)}

	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(color.Black),
		Face: basicfont.Face7x13,
		//Face: inconsolata.Regular8x16,
		Dot: point,
	}
	d.DrawString(label)
}

func drawHelper(img *image.RGBA, node *Tree, heightAbove float64, pixelsPerBranchLength float64, pixelsPerNode int, nodesPrinted *int) int {
	xStart := int(heightAbove * pixelsPerBranchLength)
	xEnd := int((heightAbove + node.BranchLength) * pixelsPerBranchLength)
	x := xEnd

	if node.Left != nil {
		yStart := drawHelper(img, node.Left, heightAbove+node.BranchLength, pixelsPerBranchLength, pixelsPerNode, nodesPrinted)
		yEnd := *nodesPrinted*pixelsPerNode + pixelsPerNode/2
		sketch.VLine(img, x, yStart, yEnd, color.Black)
	}

	y := *nodesPrinted*pixelsPerNode + pixelsPerNode/2
	sketch.HLine(img, xStart, xEnd, y, color.Black)
	*nodesPrinted = *nodesPrinted + 1
	addLabel(img, node.Name, x+5, *nodesPrinted*pixelsPerNode)

	if node.Right != nil {
		yStart := y
		yEnd := drawHelper(img, node.Right, heightAbove+node.BranchLength, pixelsPerBranchLength, pixelsPerNode, nodesPrinted)
		x := xEnd
		sketch.VLine(img, x, yStart, yEnd, color.Black)
	}

	return y
}

func Draw(node *Tree, imageWidth int, imageHeight int) (*image.RGBA, error) {
	//totalHeight := Height(node)
	img := image.NewRGBA(image.Rect(0, 0, imageWidth, imageHeight))
	sketch.FilledRectangle(img, 0, 0, imageWidth, imageHeight, color.White)

	if node == nil {
		return img, fmt.Errorf("Error: unable to draw an empty tree\n")
	} else {
		var nodesPrinted int = 0
		drawHelper(img, node, 0, 100, 10, &nodesPrinted)
		return img, nil
	}
}
