#!/bin/bash
# Cross-platform initialization script for devcontainer
# Replaces the Node.js initializeCommand for systems without Node

# Create ~/.claude directory if it doesn't exist
CLAUDE_DIR="$HOME/.claude"
if [ ! -d "$CLAUDE_DIR" ]; then
    mkdir -p "$CLAUDE_DIR"
    echo "Created $CLAUDE_DIR"
fi

# Get project name from current working directory
PROJECT_NAME=$(basename "$(pwd)")

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create .env file with COMPOSE_PROJECT_NAME
ENV_FILE="$SCRIPT_DIR/.env"
echo "COMPOSE_PROJECT_NAME=${PROJECT_NAME}-arm64" > "$ENV_FILE"

echo "Created .env with COMPOSE_PROJECT_NAME=${PROJECT_NAME}-arm64 at $ENV_FILE"
