// Package ontology provides functions for gene ontology enrichment analysis.
package ontology

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
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

	for i := range records {
		answer[records[i].Id] = &Ontology{
			Name: records[i].Name,
			Id:   records[i].Name,
		}
	}
	for i := range records {
		for j = range records[i].Parents {
			answer[records[i].Id].Parents = append(answer[records[i].Id].Parents, answer[records[i].Parents[j].Id])
		}
		for j = range records[i].Children {
			answer[records[i].Id].Children = append(answer[records[i].Id].Children, answer[records[i].Children[j].Id])
		}
	}
	return answer
}

func geneAssignmentsFromGaf(records []gaf.Gaf, terms map[string]*Ontology) {
	for i := range records {
		terms[records[i].GoId].Genes = append(terms[records[i].GoId].Genes, records[i].DbObjectSymbol)
	}
}

// genesToOntologies converts a map of GoIds to pointers to Ontology structs to
// a map[string][]*Ontology, where string is a gene or protein name and the []*Ontology is
// the set of gene-associated ontological terms.
func genesToOntologies(terms map[string]*Ontology) map[string][]*Ontology {
	var j int
	var foundInMap bool
	var answer = make(map[string][]*Ontology)

	for i := range terms {
		for j = range terms[i].Genes {
			if _, foundInMap = answer[terms[i].Genes[j]]; !foundInMap {
				answer[terms[i].Genes[j]] = []*Ontology{terms[i]}
			} else {
				answer[terms[i].Genes[j]] = append(answer[terms[i].Genes[j]], terms[i])
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
	var foundInMap bool

	for i := range filledSpace {
		currLen = filledSpace[i].ChromEnd - filledSpace[i].ChromStart
		currGene = filledSpace[i].Name
		if _, foundInMap = answerCounts[currGene]; !foundInMap {
			answerCounts[currGene] = 0
		}
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
			answer[i] += geneProportions[ontologies[i].Genes[j]]
		}
	}
	return answer
}

func ThreeDGreat(queries []bed.Bed, chromSizes map[string]chromInfo.ChromInfo, genes map[string]*gtf.Gene, contacts []bedpe.BedPe, annotations []gaf.Gaf, oboMap map[string]*obo.Obo) {
	var filledSpaceIntervals []interval.Interval
	var queryOverlaps []interval.Interval
	var kCache = make(map[string]int) //mapping Go Term Names to number of query elements that overlap. This is 'k' in the binomial test.
	var n = len(queries)              // this stores the number of query elements, which is the number of trials in the binomial test.
	var currOverlap int
	var currOverlapGene string
	var currOntology *Ontology
	var foundInMap bool

	tssBed := gtf.GenesToTssBed(genes, chromSizes, true)
	filledSpace := Fill3dSpace(contacts, tssBed, chromSizes)
	ontologies := OboToOntology(oboMap)
	geneAssignmentsFromGaf(annotations, ontologies)
	geneOntologies := genesToOntologies(ontologies)

	proportionsForGenes := geneProportionOfGenome(filledSpace)
	proportionsForTerms := termProportionOfGenome(ontologies, proportionsForGenes) // this stores the proportion of the genome that is covered by each term. Values are the 'p', or success probability, in the binomial test

	for i := range filledSpace {
		filledSpaceIntervals = append(filledSpaceIntervals, filledSpace[i])
	}
	tree := interval.BuildTree(filledSpaceIntervals)

	bed.AllToMidpoint(queries)

	for i := range queries {
		queryOverlaps = interval.Query(tree, queries[i], "any")
		for currOverlap = range queryOverlaps {
			currOverlapGene = queryOverlaps[currOverlap].(bed.Bed).Name // type assert as bed and extract name of gene assigned to query
			for currOntology = range geneOntologies[currOverlapGene] {
				if _, foundInMap = kCache[currOntology.Name]; !foundInMap {
					kCache[currOntology.Name] = 0
				}
				kCache[currOntology.Name] += 1
			}
		}
	}

	var enrichment float64
	for i := range proportionsForTerms {
		enrichment = numbers.BinomialRightSummation(n, kCache[i], proportionsForTerms[i], true)
		fmt.Printf("Enrichment: %e. Term: %s.\n", enrichment, i)
	}
}

// BIG TODO LIST
/*
1. Filter obsolete Obos and filter "Not" Gaf entries.
2 1/2: Potentially assigning genes to parent nodes in ontology tree.
6. Assign each query bed to gene
7. Calculate binomial(n, k, p) for each term
8. Return p value for each Ontology
*/
