// electron/scripts/after-sign.js
const path = require('path');
const { notarize } = require('@electron/notarize');

module.exports = async function (context) {
  // Only run for macOS builds
  if (process.platform !== 'darwin') {
    return;
  }

  if (process.env.SKIP_NOTARIZE === 'true') {
    console.log('[after-sign] SKIP_NOTARIZE=true -> skipping notarization');
    return;
  }

  const { appOutDir, packager } = context;
  const appName = packager.appInfo.productFilename; // e.g. "idepGolem"
  const appPath = path.join(appOutDir, `${appName}.app`);

  console.log(`[after-sign] Notarizing ${appPath}`);

  await notarize({
    appBundleId: 'com.idepgolem.app',
    appPath,
    appleId: process.env.APPLE_ID,
    appleIdPassword: process.env.APPLE_APP_SPECIFIC_PASSWORD,
    teamId: process.env.APPLE_TEAM_ID,
  });

  console.log('[after-sign] Notarization complete');
};
