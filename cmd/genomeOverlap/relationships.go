package main

import "fmt"

func printRelationships() {
	o := fmt.Sprintf("\t*------*\n\t    *------*\n\n\n")
	oi := fmt.Sprintf("\t    *------*\n\t*------*\n\n\n")
	d :=  fmt.Sprintf("\t  *--*\n\t*------*\n\n\n")
	di :=  fmt.Sprintf("\t*------*\n\t  *--*\n\n\n")
	m :=  fmt.Sprintf("\t*------*\n\t       *------*\n\n\n")
	mi :=  fmt.Sprintf("\t       *------*\n\t*------*\n\n\n")
	s :=  fmt.Sprintf("\t*---*\n\t*------*\n\n\n")
	si :=  fmt.Sprintf("\t*------*\n\t*---*\n\n\n")
	f :=  fmt.Sprintf("\t   *---*\n\t*------*\n\n\n")
	fi :=  fmt.Sprintf("\t*------*\n\t   *---*\n\n\n")
	gt :=  fmt.Sprintf("\t*---*\n\t        *---*\n\n\n")
	lt :=  fmt.Sprintf("\t        *---*\n\t*---*\n\n\n")
	e :=  fmt.Sprintf("\t*------*\n\t*------*\n\n\n")


	fmt.Println("Valid relationships are as follows")
	fmt.Println("Top line:    target")
	fmt.Println("Bottom line: query\n")
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
