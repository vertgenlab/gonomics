package openai

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"reflect"
	"testing"
)

func TestGetChatResponse(t *testing.T) {
	// Create a mock API server for testing
	mockServer := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// Check the request headers and body
		if r.Header.Get("Authorization") != "Bearer your-api-token" {
			t.Errorf("Error: Unexpected Authorization header")
		}
		if r.Header.Get("Content-Type") != "application/json" {
			t.Errorf("Error: Unexpected Content-Type header")
		}

		// Simulate a successful API response
		response := ChatGPTResponse{
			Choices: []string{"Hello", "Hi there"},
		}
		jsonResponse, err := json.Marshal(response)
		if err != nil {
			t.Errorf("Failed to marshal JSON response")
		}

		// Return the JSON response
		w.Header().Set("Content-Type", "application/json")
		w.WriteHeader(http.StatusOK)
		w.Write(jsonResponse)
	}))
	defer mockServer.Close()

	// Initialize the ChatGPTClient with the mock API server URL and token
	client := NewClient(mockServer.URL, "your-api-token")

	// Call the GetChatResponse function
	prompt := "Hello"
	choices, err := client.GetChatResponse(prompt)
	if err != nil {
		t.Errorf("Error: Failed to get chat response: %v", err)
	}

	// Validate the response
	expectedChoices := []string{"Hello", "Hi there"}
	if !reflect.DeepEqual(choices, expectedChoices) {
		t.Errorf("Unexpected chat response. Expected: %v, Got: %v", expectedChoices, choices)
	}
}
