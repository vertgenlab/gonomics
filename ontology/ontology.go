// Package ontology provides functions for gene ontology enrichment analysis.
package ontology

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
	"strings"
)

type Ontology struct {
	Name     string
	Id       string
	Parents  []*Ontology
	Children []*Ontology
	Genes    []string
}

// OboToOntology converts an input map[string]*obo.Obo to a map[string]*Ontology.
// returns map[string]*Ontology where string is Id.
func OboToOntology(records map[string]*obo.Obo) map[string]*Ontology {
	var j int
	var answer = make(map[string]*Ontology)
	var currOntology *Ontology

	for i := range records {
		answer[records[i].Id] = &Ontology{
			Name: records[i].Name,
			Id:   records[i].Id,
		}
	}
	for i := range records {
		currOntology = answer[records[i].Id]
		for j = range records[i].Parents {
			currOntology.Parents = append(currOntology.Parents, answer[records[i].Parents[j].Id])
		}
		for j = range records[i].Children {
			currOntology.Children = append(currOntology.Children, answer[records[i].Children[j].Id])
		}
	}
	return answer
}

// GeneAssignmentsFromGaf translates a slice of gaf structs into a map where the string key is a GO term id and the value is an ontology struct
func GeneAssignmentsFromGaf(records []gaf.Gaf, terms map[string]*Ontology) {
	records = gaf.RemoveDuplicates(records)
	for i := range records {
		if terms[records[i].GoId] != nil {
			terms[records[i].GoId].Genes = append(terms[records[i].GoId].Genes, records[i].DbObjectSymbol)
		}
		// TODO: what if records[i].GoId not found in map
	}
}

// GenesToOntologies converts a map of GoIds to pointers to Ontology structs to
// a map[string][]*Ontology, where string is a gene or protein name and the []*Ontology is
// the set of gene-associated ontological terms.
func GenesToOntologies(terms map[string]*Ontology) map[string][]*Ontology {
	var j int
	var foundInMap bool
	var answer = make(map[string][]*Ontology)
	var currGeneName string

	for i := range terms {
		for j = range terms[i].Genes {
			currGeneName = terms[i].Genes[j]
			if _, foundInMap = answer[currGeneName]; !foundInMap {
				answer[currGeneName] = []*Ontology{terms[i]}
			} else {
				answer[currGeneName] = append(answer[currGeneName], terms[i])
			}
		}
	}

	return answer
}

// geneProportionOfGenome returns the proportion of genome occupied by each gene.
// answer is of form map[string]float64, where the key string is the gene name and the value is the proportion
func geneProportionOfGenome(filledSpace []bed.Bed) map[string]float64 {
	var answerCounts = make(map[string]int)
	var answer = make(map[string]float64)
	var totalBases, currLen int = 0, 0
	var currGene string

	for i := range filledSpace {
		currLen = filledSpace[i].ChromEnd - filledSpace[i].ChromStart
		currGene = strings.ToUpper(filledSpace[i].Name)
		answerCounts[currGene] += currLen
		totalBases += currLen
	}

	for i := range answerCounts {
		answer[i] = float64(answerCounts[i]) / float64(totalBases)
	}

	return answer
}

// termProportionOfGenome returns a map[string]float64, where keys are the names of a gene ontology term and values are the
// proportion of the genome covered by that term.
func termProportionOfGenome(ontologies map[string]*Ontology, geneProportions map[string]float64) map[string]float64 {
	var answer = make(map[string]float64)
	var j int
	for i := range ontologies {
		answer[i] = 0
		for j = range ontologies[i].Genes {
			answer[i] += geneProportions[strings.ToUpper(ontologies[i].Genes[j])]
		}
	}
	return answer
}

// ThreeDGreat uses all inputs to determine which areas of the genome, based on 3D genome contacts, are associated with
// a gene, and therefore it's function, and will output the enrichment score for a given GO term.
// This function takes in a query file as a bed struct, a chrom.sizes file, a map of genes in gtf form, a contact file in
// bedpe form, an annotation file in gaf form and a map that takes GO term ID's and relates them to other GO term features
// in obo format, and finally it takes a string which will write out all regions of the genome with their assigned
// genes in the name column and their assigned ontologies in the annotation field.
func ThreeDGreat(queries []bed.Bed, chromSizes map[string]chromInfo.ChromInfo, geneFile string, contacts []bedpe.BedPe, annotations []gaf.Gaf, oboMap map[string]*obo.Obo, out3dOntology string, geneEnrichments bool, termEnrichments bool) {
	var err error
	var filledSpace, tssBed []bed.Bed
	var filledSpaceIntervals []interval.Interval
	var queryOverlaps []interval.Interval
	var kCache = make(map[string]int) //mapping Go Term Names to number of query elements that overlap. This is 'k' in the binomial test.
	var n = len(queries)              // this stores the number of query elements, which is the number of trials in the binomial test.
	var currOverlapGene string
	var currOntologyIndex int
	var currOntologyName string
	name := strings.TrimSuffix(out3dOntology, ".bed")
	var ontologiesForCurrGene []*Ontology

	geneString := strings.Split(geneFile, ".")
	if geneString[len(geneString)-1] != "bed" {
		genes := gtf.Read(geneFile)
		tssBed = gtf.GenesToTssBed(genes, chromSizes, true)
	} else {
		tssBed = bed.Read(geneFile)
	}
	bed.SortByCoord(tssBed)
	filledSpace = Fill3dSpace(contacts, tssBed, chromSizes)
	ontologies := OboToOntology(oboMap)
	GeneAssignmentsFromGaf(annotations, ontologies)
	geneOntologies := GenesToOntologies(ontologies)

	//TODO: filtering out genes without ontology annotations? as an option?

	if out3dOntology != "" {
		write3dOntologies(out3dOntology, geneOntologies, filledSpace)
	}

	var proportionsForGenes map[string]float64
	proportionsForGenes = geneProportionOfGenome(filledSpace)
	if geneEnrichments {
		geneOut := fileio.EasyCreate(name + ".geneProportions.txt")
		_, err = fmt.Fprintf(geneOut, "Gene\tProportion\n")
		exception.PanicOnErr(err)
		for i := range proportionsForGenes {
			_, err = fmt.Fprintf(geneOut, "%s\t%e\n", i, proportionsForGenes[i])
			exception.PanicOnErr(err)
		}
		err = geneOut.Close()
		exception.PanicOnErr(err)
	}

	for i := range filledSpace {
		filledSpaceIntervals = append(filledSpaceIntervals, filledSpace[i])
	}
	tree := interval.BuildTree(filledSpaceIntervals)

	bed.AllToMidpoint(queries)

	for i := range queries {
		queryOverlaps = interval.Query(tree, queries[i], "any")
		if len(queryOverlaps) != 1 {
			log.Fatalf("Query overlapped multiple regions in filled space. %v", queries[i])
		}
		currOverlapGene = queryOverlaps[0].(bed.Bed).Name // type assert as bed and extract name of gene assigned to query
		ontologiesForCurrGene = geneOntologies[currOverlapGene]
		for currOntologyIndex = range ontologiesForCurrGene {
			currOntologyName = ontologiesForCurrGene[currOntologyIndex].Id
			kCache[currOntologyName] += 1
		}
	}

	var proportionsForTerms map[string]float64
	var enrichment float64
	if termEnrichments {
		proportionsForTerms = termProportionOfGenome(ontologies, proportionsForGenes) // this stores the proportion of the genome that is covered by each term. Values are the 'p', or success probability, in the binomial test
		out := fileio.EasyCreate(name + ".termProportions.txt")
		_, err = fmt.Fprintf(out, "Term\tName\tProportion\n")
		enrichOut := fileio.EasyCreate(name + ".termEnrichment.txt")
		_, err = fmt.Fprintf(enrichOut, "Term\tName\tEnrichment\n")
		for i := range proportionsForTerms {
			if proportionsForTerms[i] > 0 {
				_, err = fmt.Fprintf(out, "%s\t%s\t%e\n", i, ontologies[i].Name, proportionsForTerms[i])
				exception.PanicOnErr(err)
				enrichment = numbers.BinomialRightSummation(n, kCache[i], proportionsForTerms[i], true)
				_, err = fmt.Fprintf(enrichOut, "%s\t%s\t%e\n", i, ontologies[i].Name, enrichment)
				exception.PanicOnErr(err)
			} else {
				continue
			}
		}
		err = out.Close()
		exception.PanicOnErr(err)
		err = enrichOut.Close()
		exception.PanicOnErr(err)
	}
}

// write3dOntologies take a 3D filled space bed, a map of gene names to their ontologies and an output file name and
// write out a bed that contains all the filled space information as well as the ontologies that relate to the gene assignments in teh Annotation field of bed
func write3dOntologies(filename string, geneToOnt map[string][]*Ontology, filledSpace []bed.Bed) {
	var onts []string

	for i := range filledSpace {
		filledSpace[i].FieldsInitialized = 8
		onts = ontologiesToStrings(geneToOnt[filledSpace[i].Name])
		filledSpace[i].FieldsInitialized += len(onts)
		filledSpace[i].Score = filledSpace[i].Score + 0 //just to make sure it isn't empty, or overwritten
		filledSpace[i].Strand = '.'
		filledSpace[i].Annotation = append(filledSpace[i].Annotation, onts...)
	}
	bed.Write(filename, filledSpace)
}

// ontologiesToStrings take a slice of pointers to Ontology structs and stores their name field in a slice of strings
func ontologiesToStrings(onts []*Ontology) []string {
	var answer []string

	for o := range onts {
		answer = append(answer, onts[o].Name)
	}

	return answer
}

// BIG TODO LIST
/*
1. Filter obsolete Obos and filter "Not" Gaf entries.
2. Potentially assigning genes to parent nodes in ontology tree.

Note: gene slice must have unique entries.

*/
