package apis

import (
	"context"
	"errors"
	"fmt"
	"log"
	"os"

	"github.com/google/generative-ai-go/genai"
	"google.golang.org/api/option"
)

func NewClient(ctx context.Context) *genai.Client{
    apiKey, ok := os.LookupEnv("API_KEY")
    if !ok || apiKey == "" {
        log.Panic( errors.New("API_KEY environment variable not found or empty"))
    }

    client, err := genai.NewClient(ctx, option.WithAPIKey(apiKey))
    if err != nil {
        log.Panic(fmt.Errorf("failed to create genai client: %v", err))
    }

    return client
}