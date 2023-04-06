// Package sketch provides utilities and color palettes for plotting graphs and images
package sketch

import (
	"image"
	"image/color"

	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"
	// "golang.org/x/image/font/inconsolata"
)

func HLine(img *image.RGBA, xStart int, xEnd int, y int, col color.Color) {
	for x := xStart; x < xEnd; x++ {
		img.Set(x, y, col)
	}
}

func VLine(img *image.RGBA, x int, yStart int, yEnd int, col color.Color) {
	for y := yStart; y < yEnd; y++ {
		img.Set(x, y, col)
	}
}

func Rectangle(img *image.RGBA, xOne int, yOne int, xTwo int, yTwo int, col color.Color) {
	HLine(img, xOne, xTwo, yOne, col)
	HLine(img, xOne, xTwo, yTwo, col)
	VLine(img, xOne, yOne, yTwo, col)
	VLine(img, xTwo, yOne, yTwo, col)
}

func FilledRectangle(img *image.RGBA, xOne int, yOne int, xTwo int, yTwo int, col color.Color) {
	for x := xOne; x < xTwo; x++ {
		for y := yOne; y < yTwo; y++ {
			img.Set(x, y, col)
		}
	}
}

func Text(img *image.RGBA, label string, xStart int, yStart int) {
	point := fixed.Point26_6{fixed.Int26_6(xStart * 64), fixed.Int26_6(yStart * 64)}

	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(color.Black),
		Face: basicfont.Face7x13,
		// Face: inconsolata.Regular8x16,
		Dot: point,
	}
	d.DrawString(label)
}
