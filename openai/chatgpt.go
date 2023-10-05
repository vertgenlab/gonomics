// Package openai provides a client for interfacing with the OpenAI ChatGPT service.
package openai

import (
	"bytes"
	"encoding/json"
	"log"
	"net/http"
	"os"
)

// ChatGPTClient is the primary struct representing a client to interact with OpenAI's ChatGPT.
type ChatGPTClient struct {
	APIURL string       // The API endpoint URL
	Token  string       // API token for authentication
	Client *http.Client // HTTP client used for sending requests
}

// ChatGPTResponse represents the response structure returned by the OpenAI ChatGPT API.
type ChatGPTResponse struct {
	Choices []struct {
		Message struct {
			Role    string `json:"role"`    // Role of the message sender (either "system", "user", or "assistant")
			Content string `json:"content"` // The content of the message
		} `json:"message"`
	} `json:"choices"`
}

// NewClient creates and returns a new ChatGPTClient given the API URL and authentication token.
func NewClient(apiURL, token string) *ChatGPTClient {
	return &ChatGPTClient{
		APIURL: apiURL,
		Token:  token,
		Client: &http.Client{},
	}
}

// GetAPIToken retrieves the API token from the environment using the provided apiKey variable name.
func GetAPIToken(apiKey string) string {
	return os.Getenv(apiKey)
}

// PostChatResponse sends a series of messages to the OpenAI ChatGPT API and returns the assistant's response.
// It takes a slice of messages and constructs the API request. If there's an error at any step, it terminates the program.
func (c *ChatGPTClient) PostChatResponse(messages []interface{}) string {
	jsonBody, err := json.Marshal(map[string]interface{}{
		"model":      gptThreeFiveTurbo,
		"messages":   messages,
		"max_tokens": 248,
	})
	if err != nil {
		log.Fatalf("Error: Encoding request body: %v", err)
	}

	req, err := http.NewRequest("POST", c.APIURL, bytes.NewBuffer(jsonBody))
	if err != nil {
		log.Fatalf("Error: Creating the request: %v", err)
	}
	req.Header.Set("Content-Type", "application/json")
	req.Header.Set("Authorization", "Bearer "+c.Token)

	resp, err := c.Client.Do(req)
	if err != nil {
		log.Fatalf("Error: Sending the request: %v", err)
	}
	defer resp.Body.Close()

	var response ChatGPTResponse
	err = json.NewDecoder(resp.Body).Decode(&response)
	if err != nil {
		log.Fatalf("Error: Decoding JSON response: %v", err)
	}

	return response.Choices[0].Message.Content
}
