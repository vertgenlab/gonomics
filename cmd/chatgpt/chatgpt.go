package main

import (
	"bufio"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/openai"
)

func main() {
	apiKey := openai.GetAPIToken("OPENAI_API_TOKEN")
	if apiKey == "" {
		log.Fatalf("Error: API Token OPENAI_API_KEY=\"\" not set")
	}

	client := openai.NewClient(openai.ApiEndpoint, apiKey)

	messages := []interface{}{
		map[string]interface{}{
			"role":    "system",
			"content": "You are a helpful assistant.",
		},
	}

	reader := bufio.NewReader(os.Stdin)

	for {
		fmt.Print("You: ")
		userInput, err := reader.ReadString('\n')
		if err != nil {
			log.Fatalf("Error reading input: %v", err)
		}

		// Remove the trailing newline
		userInput = userInput[:len(userInput)-1]

		if userInput == "exit" {
			break
		}

		// Append the user's message
		messages = append(messages, map[string]interface{}{"role": "user", "content": userInput})

		// Get the assistant's response and then append it
		response := client.PostChatResponse(messages)
		messages = append(messages, map[string]interface{}{"role": "assistant", "content": response})

		fmt.Println("Assistant:", response)
	}
}
