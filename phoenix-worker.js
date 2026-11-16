/* phoenix-worker.js — runs analyses off the main thread.
 *
 * The engine modules are plain UMD scripts, so a classic worker can load them
 * with importScripts and call exactly the same code the page uses. Results are
 * structured-cloneable (plain objects, strings, numbers); interactive sessions
 * hold functions and are therefore run on the main thread instead.
 */
importScripts(
  'phoenix-core.js',
  'phoenix-models.js',
  'phoenix-linalg.js',
  'phoenix-infer.js',
  'phoenix-support.js',
  'phoenix-continuous.js',
  'phoenix-interactive.js',
  'phoenix-programs.js'
);

self.onmessage = function (event) {
  const { id, programId, text, opts, trees } = event.data;
  try {
    const programs = self.PPrograms;
    const core = self.PCore;
    const program = programs.get(programId);
    const ctx = { text: text || '', opts: opts || {} };
    if (trees && trees.length) ctx.trees = trees.map(t => core.parseNewick(t));
    const result = program.run(ctx);
    if (result.session) throw new Error('interactive programs must run on the main thread');
    self.postMessage({ id, ok: true, result });
  } catch (error) {
    self.postMessage({ id, ok: false, error: error && error.message ? error.message : String(error) });
  }
};
