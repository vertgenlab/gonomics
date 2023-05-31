package openai

import (
	"bytes"
	"encoding/json"
	"net/http"
	"os"
)

type ChatGPTClient struct {
	APIURL string
	Token  string
}

type ChatGPTResponse struct {
	Choices []string `json:"choices"`
}

func NewClient(apiURL, token string) *ChatGPTClient {
	return &ChatGPTClient{
		APIURL: apiURL,
		Token:  token,
	}
}

// GetAPIToken retrieves the API token from an environment variable.
// It returns the token if found, or an empty string if not found.
func GetAPIToken(apiKey string) string {
	token := os.Getenv(apiKey)
	return token
}

func (c *ChatGPTClient) GetChatResponse(prompt string) ([]string, error) {
	url := c.APIURL + "/chat/completions"
	requestBody := map[string]interface{}{
		"prompt": prompt,
	}
	jsonBody, err := json.Marshal(requestBody)
	if err != nil {
		return nil, err
	}

	req, err := http.NewRequest("POST", url, bytes.NewBuffer(jsonBody))
	if err != nil {
		return nil, err
	}

	req.Header.Set("Content-Type", "application/json")
	req.Header.Set("Authorization", "Bearer "+c.Token)

	client := &http.Client{}
	resp, err := client.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	var response ChatGPTResponse
	err = json.NewDecoder(resp.Body).Decode(&response)
	if err != nil {
		return nil, err
	}

	return response.Choices, nil
}
