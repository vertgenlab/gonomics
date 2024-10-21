package main

// Example of JSON layout can be found in the link below
// http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/chr1%3A878884%3AC%3AT,chr1%3A878917%3AT%3AA,chr1%3A878920%3AT%3AA,chr1%3A878991%3AGTGTT%3AG,chr1%3A879229%3AA%3AT,chr1%3A879231%3AA%3AC,chr1%3A879897%3AT%3AC,chr1%3A879957%3AG%3AT/annotation?assembly=grch38

type Responses struct { // the response for each variant queried
	Responses []Response `json:"response"`
}

type Response struct { // all the results for a single variant (e.g. multiple transcripts)
	Results []Result `json:"result"`
}

type Result struct { // data on a particular variant on a particular transcript
	Chr             string          `json:"chromosome"`
	Start           int             `json:"start"`
	Ref             string          `json:"reference"`
	Alt             string          `json:"alternate"`
	SnpId           string          `json:"id"`
	ConsequenceType string          `json:"displayConsequenceType"`
	Consequences    []Consequence   `json:"consequenceTypes"`
	PopAlleleFreqs  []PopAlleleFreq `json:"populationFrequencies"`
	// Fields below exist but are not collected for this cmd
	// conservation
	// geneExpression
	// geneTraitAssociation
	// geneDrugInteraction
	// cytoband
	// repeat
}

type PopAlleleFreq struct { // pop allele frequency from gnomAD and 1kgp
	Study      string  `json:"study"`
	Population string  `json:"population"`
	RefAf      float64 `json:"refAlleleFreq"`
	AltAf      float64 `json:"altAlleleFreq"`
}

type Consequence struct { // information about variant effect
	GeneName          string            `json:"geneName"`
	GeneId            string            `json:"ensemblGeneId"`
	TranscriptId      string            `json:"ensemblTranscriptId"`
	Strand            string            `json:"strand"`
	Biotype           string            `json:"biotype"`
	ProteinAnnotation ProteinAnnotation `json:"proteinVariantAnnotation"`
}

type ProteinAnnotation struct { // variants effect on protein
	Pos                int                  `json:"position"`
	Ref                string               `json:"reference"`
	Alt                string               `json:"alternate"`
	SubstitutionScores []SubstitutionScores `json:"substitutionScores"`
}

type SubstitutionScores struct { // sift & polyphen scores
	Source      string  `json:"source"`
	Score       float64 `json:"score"`
	Description string  `json:"description"`
}
