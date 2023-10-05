package openai

import (
	"bytes"
	"io/ioutil"
	"net/http"
	"os"
	"testing"
)

// MockRoundTripper struct to mock HTTP calls
type MockRoundTripper struct {
	resp *http.Response
	err  error
}

// RoundTrip is the mock implementation of http.RoundTripper interface
func (m *MockRoundTripper) RoundTrip(req *http.Request) (*http.Response, error) {
	return m.resp, m.err
}

func TestPostChatResponse(t *testing.T) {
	client := NewClient("http://mockapi", "mocktoken")

	// Creating a mocked HTTP response
	body := `{
		"choices": [{
			"message": {
				"role": "assistant",
				"content": "Hello, world!"
			}
		}]
	}`
	rt := &MockRoundTripper{
		resp: &http.Response{
			StatusCode: 200,
			Body:       ioutil.NopCloser(bytes.NewReader([]byte(body))),
		},
	}
	client.Client.Transport = rt

	messages := []interface{}{
		map[string]interface{}{
			"role":    "system",
			"content": "You are a helpful assistant.",
		},
	}

	resp := client.PostChatResponse(messages)

	if resp != "Hello, world!" {
		t.Errorf("Expected response 'Hello, world!', got %s", resp)
	}
}

func TestGetAPIToken(t *testing.T) {
	const mockToken = "mocktokenvalue"
	_ = os.Setenv("MOCK_API_TOKEN", mockToken) // Ignoring error for simplicity
	defer os.Unsetenv("MOCK_API_TOKEN")

	token := GetAPIToken("MOCK_API_TOKEN")

	if token != mockToken {
		t.Errorf("Expected token '%s', got %s", mockToken, token)
	}
}
