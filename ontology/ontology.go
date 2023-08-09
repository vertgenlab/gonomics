// Package ontology provides functions for gene ontology enrichment analysis.
package ontology

import (
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

// BIG TODO LIST
/*
1. Filter obsolete Obos and filter "Not" Gaf entries.
2 1/2: Potentially assigning genes to parent nodes in ontology tree.
5. Calculate # of bases for each ontology (this is the binomial parameter p)
6. Assign each query bed to gene
7. Calculate binomial(n, k, p) for each term
8. Return p value for each Ontology
*/
