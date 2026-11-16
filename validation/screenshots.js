#!/usr/bin/env node
/* validation/screenshots.js — capture screenshots of the running application.
 * Requires Playwright. Writes PNGs into docs/screenshots/.
 *   npm run screenshots
 */
'use strict';
const fs = require('fs');
const path = require('path');
let chromium;
try { ({ chromium } = require('playwright')); }
catch (_) {
  try { ({ chromium } = require('/home/user/.npm-global/lib/node_modules/playwright')); }
  catch (__) { console.log('SKIP: Playwright not installed.'); process.exit(0); }
}

const ROOT = path.join(__dirname, '..');
const OUT = path.join(ROOT, 'docs', 'screenshots');
const URL = 'file://' + path.join(ROOT, 'index.html');

(async () => {
  fs.mkdirSync(OUT, { recursive: true });
  const browser = await chromium.launch();

  async function shot(page, name, fullPage) {
    await page.screenshot({ path: path.join(OUT, name + '.png'), fullPage: !!fullPage });
    console.log('wrote docs/screenshots/' + name + '.png');
  }

  const page = await browser.newPage({ viewport: { width: 1440, height: 950 } });
  await page.goto(URL);
  await page.waitForTimeout(500);
  await shot(page, '01-landing');

  await page.locator('#demos').scrollIntoViewIfNeeded();
  await page.waitForTimeout(300);
  await shot(page, '02-demos');

  await page.click('.demo-run[data-demo="bootstrap"]');
  await page.waitForTimeout(5000);
  await page.locator('#suite-output').scrollIntoViewIfNeeded();
  await page.waitForTimeout(400);
  await shot(page, '03-bootstrap-consensus');

  await page.click('.program-item[data-program="MOVE"]');
  await page.waitForTimeout(200);
  await page.click('#runProgram');
  await page.waitForTimeout(1200);
  await shot(page, '04-interactive-move');

  await page.click('.program-item[data-program="DNADIST"]');
  await page.waitForTimeout(100);
  await page.click('#runProgram');
  await page.waitForTimeout(800);
  await shot(page, '05-distance-matrix');

  const mobile = await browser.newPage({ viewport: { width: 390, height: 844 } });
  await mobile.goto(URL);
  await mobile.waitForTimeout(500);
  await shot(mobile, '06-mobile-landing');

  const tablet = await browser.newPage({ viewport: { width: 820, height: 1180 } });
  await tablet.goto(URL);
  await tablet.waitForTimeout(500);
  await shot(tablet, '07-tablet');

  await browser.close();
})();
