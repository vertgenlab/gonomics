package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/openai"
)

func main() {
	// Define command-line flags
	tokenFlag := flag.String("token", "", "API token")
	questionFlag := flag.String("question", "", "Question to ask the Chatbot")

	flag.Parse()
	var apiToken string

	// Check if API token is provided
	if *tokenFlag == "" {
		// Attempt to retrieve API token from environment variable
		apiToken := os.Getenv("OPENAI_API_TOKEN")
		if apiToken != "" {
			tokenFlag = &apiToken
		} else {
			log.Fatal(
				"Error: API token is missing. Please provide it using the -token flag or set GPT_API_TOKEN environment variable.",
			)
		}
	}
	*tokenFlag = apiToken

	// Create ChatGPT client
	client := openai.NewClient("https://api.openai.com/v1", *tokenFlag)

	// Check if question is provided
	if *questionFlag == "" {
		log.Fatal("Error: Question is missing. Please provide it using the -question flag.")
	}

	// Generate a response
	responses, err := client.GetChatResponse(*questionFlag)
	if err != nil {
		log.Fatal("Failed to generate a response:", err)
	}

	fmt.Println("Response:")
	for _, response := range responses {
		fmt.Println(response)
	}
}
