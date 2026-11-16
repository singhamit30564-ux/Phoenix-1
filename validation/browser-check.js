#!/usr/bin/env node
/* validation/browser-check.js — headless browser checks for the Phoenix UI.
 *
 * Requires Playwright (not bundled with the app; the app itself has no
 * dependencies). If Playwright is unavailable the script exits 0 with a clear
 * "skipped" message rather than failing the build.
 *
 *   npm run browser-check
 *
 * Writes validation/browser-report.json.
 */
'use strict';

const fs = require('fs');
const path = require('path');

let chromium;
try {
  ({ chromium } = require('playwright'));
} catch (_) {
  try {
    ({ chromium } = require('/home/user/.npm-global/lib/node_modules/playwright'));
  } catch (__) {
    console.log('SKIP: Playwright is not installed; browser checks were not run.');
    console.log('      Install with: npm i -g playwright && npx playwright install chromium');
    process.exit(0);
  }
}

const ROOT = path.join(__dirname, '..');
const URL = 'file://' + path.join(ROOT, 'index.html');

const SAMPLE_TREES = {
  CONTRAST: { program: 'numbers', tree: '((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3,E:0.5);' },
  RETREE: { program: 'discrete', tree: '((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);', outgroup: 'D' },
  CONSENSE: { program: 'discrete', tree: null },
  TREEDIST: { program: 'discrete', tree: null }
};

const results = [];

function record(name, pass, detail) {
  results.push({ name, pass, detail: detail || '' });
}

async function overflow(page) {
  return page.evaluate(() => ({
    docWidth: document.documentElement.scrollWidth,
    winWidth: window.innerWidth
  }));
}

(async () => {
  const browser = await chromium.launch();
  const context = await browser.newContext({ viewport: { width: 1280, height: 900 } });
  const page = await context.newPage();
  const consoleErrors = [];
  const pageErrors = [];
  const failedRequests = [];
  page.on('console', m => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', e => pageErrors.push(e.message));
  page.on('requestfailed', r => failedRequests.push(r.url()));

  await page.goto(URL);
  await page.waitForTimeout(400);

  // ---- 1. every program runs and renders -------------------------------
  const ids = await page.$$eval('.program-item', els => els.map(e => e.dataset.program));
  record('catalogue lists 35 programs', ids.length === 35, `found ${ids.length}`);

  let programFailures = [];
  for (const id of ids) {
    await page.click(`.program-item[data-program="${id}"]`);
    await page.waitForTimeout(60);
    const sample = SAMPLE_TREES[id];
    if (id === 'CONTRAST') {
      await page.fill('#programTree', sample.tree);
      await page.fill('#programInput', '5 2\nA 0.5 1.2\nB 0.6 1.1\nC 0.2 0.4\nD 0.25 0.45\nE 0.6 0.9');
    }
    if (id === 'RETREE') {
      await page.fill('#programTree', sample.tree);
      await page.fill('#programOptions input[data-option="outgroup"]', 'D');
    }
    if (id === 'DNAINVAR') {
      await page.fill('#programInput', '>A\nACGTACGTACGTACGT\n>C\nACGTACGTACGTACGT\n>B\nACGTACGTACGTACGA\n>D\nACGTACGTACGTACGA');
    }
    await page.click('#runProgram');
    await page.waitForTimeout(id.includes('ML') || id === 'SEQBOOT' ? 2500 : 700);
    const error = (await page.locator('#programError').textContent()).trim();
    const hasView = (await page.locator('#programView').innerHTML()).length > 60;
    const hasText = (await page.locator('#programText').textContent()).trim().length > 20;
    const metaVisible = await page.locator('#resultMeta').isVisible();
    if (error) programFailures.push(`${id}: ${error}`);
    else if (!hasView && !hasText) programFailures.push(`${id}: no output`);
    else if (!metaVisible) programFailures.push(`${id}: no metadata strip`);
  }
  record('all 35 programs render a result with metadata', programFailures.length === 0, programFailures.join('; '));

  // ---- 2. demos ---------------------------------------------------------
  const demoFailures = [];
  for (const demo of ['dna', 'parsimony', 'bootstrap']) {
    await page.click(`.demo-run[data-demo="${demo}"]`);
    await page.waitForTimeout(demo === 'bootstrap' ? 5000 : 1500);
    const svg = await page.locator('#programView svg').count();
    const text = (await page.locator('#programText').textContent()).trim();
    const err = (await page.locator('#programError').textContent()).trim();
    if (err || svg === 0 || text.length < 20) demoFailures.push(`${demo}: err=${err} svg=${svg}`);
  }
  record('three guided demos run end to end', demoFailures.length === 0, demoFailures.join('; '));

  // ---- 3. exports -------------------------------------------------------
  await page.click('.program-item[data-program="DNADIST"]');
  await page.waitForTimeout(80);
  await page.click('#runProgram');
  await page.waitForTimeout(600);
  const csvEnabled = await page.locator('#downloadCsv').isEnabled();
  const jsonEnabled = await page.locator('#downloadJson').isEnabled();
  record('distance run enables CSV and JSON export', csvEnabled && jsonEnabled, `csv=${csvEnabled} json=${jsonEnabled}`);

  // ---- 4. repeated execution resets state -------------------------------
  const firstRun = await page.locator('#programText').textContent();
  await page.click('#runProgram');
  await page.waitForTimeout(600);
  const secondRun = await page.locator('#programText').textContent();
  record('running the same program twice is stable', firstRun === secondRun);

  // ---- 5. error handling -------------------------------------------------
  await page.click('.program-item[data-program="DNAML"]');
  await page.waitForTimeout(80);
  await page.fill('#programInput', 'this is not a sequence file');
  await page.click('#runProgram');
  await page.waitForTimeout(700);
  const errorText = (await page.locator('#programError').textContent()).trim();
  record('malformed input produces a readable error', errorText.length > 10, errorText);

  // ---- 6. responsive layouts ---------------------------------------------
  await page.click('.program-item[data-program="NEIGHBOR"]');
  await page.waitForTimeout(80);
  await page.click('#runProgram');
  await page.waitForTimeout(800);
  const viewports = [
    { name: 'mobile', width: 390, height: 844 },
    { name: 'tablet', width: 820, height: 1180 },
    { name: 'desktop', width: 1440, height: 900 }
  ];
  const responsive = [];
  for (const vp of viewports) {
    await page.setViewportSize({ width: vp.width, height: vp.height });
    await page.waitForTimeout(250);
    const { docWidth, winWidth } = await overflow(page);
    const slack = docWidth - winWidth;
    responsive.push(`${vp.name}: overflow ${slack}px`);
    if (slack > 2) record(`no horizontal overflow at ${vp.name} (${vp.width}px)`, false, `overflow ${slack}px`);
    else record(`no horizontal overflow at ${vp.name} (${vp.width}px)`, true, `overflow ${slack}px`);
  }

  // ---- 7. console hygiene -------------------------------------------------
  record('no uncaught page errors', pageErrors.length === 0, pageErrors.join('; '));
  record('no console errors', consoleErrors.length === 0, consoleErrors.join('; '));
  record('no failed network requests', failedRequests.length === 0, failedRequests.join('; '));

  await browser.close();

  const passed = results.filter(r => r.pass).length;
  const report = {
    suite: 'Phoenix browser checks',
    generatedAt: new Date().toISOString(),
    url: URL,
    summary: { checks: results.length, passed, failed: results.length - passed, viewports: responsive },
    results
  };
  fs.writeFileSync(path.join(__dirname, 'browser-report.json'), JSON.stringify(report, null, 2) + '\n');

  for (const r of results) {
    console.log(`${r.pass ? 'PASS' : 'FAIL'}  ${r.name}${r.pass ? '' : '  — ' + r.detail}`);
  }
  console.log(`\n${passed}/${results.length} browser checks passed.`);
  process.exit(passed === results.length ? 0 : 1);
})().catch(error => {
  console.error('Browser check crashed:', error);
  process.exit(1);
});
