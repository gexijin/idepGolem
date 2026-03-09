#!/usr/bin/env node
/**
 * Optional: Extracts fields from host's ~/.claude.json for the container.
 *
 * The default postCreateCommand (`echo '{"hasCompletedOnboarding": true}' > ~/.claude.json`)
 * is sufficient for authentication. This script is for users who want additional
 * host config fields (e.g., oauthAccount, userID) in their container.
 *
 * To use this script:
 * 1. Add mount to devcontainer.json:
 *    "source=${localEnv:USERPROFILE}/.claude.json,target=/tmp/host-claude.json,type=bind,readonly"
 * 2. Change postCreateCommand to:
 *    "node /workspaces/vibe/.devcontainer/setup-claude-config.js /tmp/host-claude.json ~/.claude.json"
 *
 * Usage: node setup-claude-config.js <source> <target>
 *   source: path to mounted host config (e.g., /tmp/host-claude.json)
 *   target: path to write container config (e.g., ~/.claude.json)
 *
 * Note: This does NOT fix the "Claude API" UI label - that appears to be
 * determined by the authentication method, not config fields.
 */

const fs = require('fs');
const path = require('path');

const SOURCE = process.argv[2] || '/tmp/host-claude.json';
const TARGET = process.argv[3] || path.join(process.env.HOME, '.claude.json');

// Fields to extract from host config
const FIELDS_TO_COPY = [
  // Auth & account
  'oauthAccount',
  'userID',
  'hasAvailableSubscription',
  'hasOpusPlanDefault',
  'claudeCodeFirstTokenDate',
  // Onboarding state
  'hasCompletedOnboarding',
  'lastOnboardingVersion',
  'shiftEnterKeyBindingInstalled',
  'installMethod',
  'firstStartTime',
  // Misc
  'isQualifiedForDataSharing',
  'autoUpdates',
];

function main() {
  // Always set these defaults
  const containerConfig = {
    hasCompletedOnboarding: true,
    shiftEnterKeyBindingInstalled: true,
  };

  // Try to read host config
  if (fs.existsSync(SOURCE)) {
    try {
      const hostConfig = JSON.parse(fs.readFileSync(SOURCE, 'utf8'));

      // Copy specified fields if they exist
      for (const field of FIELDS_TO_COPY) {
        if (hostConfig[field] !== undefined) {
          containerConfig[field] = hostConfig[field];
        }
      }

      console.log('Extracted fields from host config:', Object.keys(containerConfig).join(', '));
    } catch (err) {
      console.warn('Warning: Could not parse host config:', err.message);
      console.log('Using defaults only');
    }
  } else {
    console.log('Host config not found at', SOURCE);
    console.log('Using defaults only');
  }

  // Write container config
  fs.writeFileSync(TARGET, JSON.stringify(containerConfig, null, 2));
  console.log('Wrote container config to', TARGET);
}

main();
