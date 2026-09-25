// Checks the actual MIT Unipept Visualizations DOM and exported figure attribution.
// https://github.com/unipept/unipept-visualizations; doi:10.1093/bioinformatics/btab590.
import { chromium } from "playwright";
import fs from "node:fs";
import path from "node:path";
if (process.argv.length !== 4)
  throw new Error(
    "Usage: node check_taxonomy.mjs REPORT/index.html NEW_OUTPUT",
  );
const out = path.resolve(process.argv[3]);
fs.mkdirSync(out, { recursive: false });
const browser = await chromium.launch({
  headless: true,
  args: ["--no-sandbox", "--disable-dev-shm-usage"],
});
const page = await browser.newPage({
  viewport: { width: 1280, height: 1100 },
  acceptDownloads: true,
});
const errors = [],
  remote = [];
page.on("pageerror", (e) => errors.push(e.message));
await page.route(/^https?:/, (route) => {
  remote.push(route.request().url());
  return route.abort();
});
await page.goto("file://" + path.resolve(process.argv[2]));
await page.waitForSelector("body[data-ready=true]", { timeout: 20000 });
for (const view of ["sunburst", "treemap", "treeview"]) {
  await page.locator(`[data-view=${view}]`).click();
  await page.waitForTimeout(200);
  await page
    .locator("#taxonomy-panel")
    .screenshot({ path: path.join(out, view + ".png") });
  if (
    (await page
      .locator(
        view === "treemap" ? "#taxonomy-stage .node" : "#taxonomy-stage > svg",
      )
      .count()) < 1
  )
    throw new Error("Missing " + view + " view");
  const measured = await page.evaluate((view) => {
    const elements = [
      ...document.querySelectorAll(
        view === "treemap" ? "#taxonomy-stage .node" : "#taxonomy-stage .arc",
      ),
    ];
    const root = elements
      .map((e) => e.__data__)
      .find((d) => d && d.data && d.data.id === 1);
    return root
      ? {
          actual: root.value,
          expected: Number(document.body.dataset.expectedTotal),
        }
      : null;
  }, view);
  if (
    measured &&
    Math.abs(measured.actual - measured.expected) >
      1e-6 * Math.max(1, measured.expected)
  )
    throw new Error("View weight differs " + JSON.stringify(measured));
  for (const metric of ["signal", "peptides"]) {
    await page.selectOption("#metric", metric);
    for (const sample of ["0", "1"]) {
      await page.selectOption("#sample", sample);
      const mass = await page.evaluate((view) => {
        const elements = [...document.querySelectorAll("#taxonomy-stage *")];
        const root = elements
          .map((e) => e.__data__)
          .find((d) => d && d.data && d.data.id === 1);
        return root
          ? {
              actual: view === "treeview" ? root.data.count : root.value,
              expected: Number(document.body.dataset.expectedTotal),
            }
          : null;
      }, view);
      if (
        !mass ||
        !Number.isFinite(mass.actual) ||
        Math.abs(mass.actual - mass.expected) >
          1e-6 * Math.max(1, mass.expected)
      )
        throw new Error(
          "View mass mismatch " +
            view +
            " " +
            sample +
            " " +
            metric +
            " " +
            JSON.stringify(mass),
        );
    }
  }
  await page.selectOption("#sample", "0");
  {
    const promise = page.waitForEvent("download");
    await page.locator("#save-svg").click();
    const download = await promise;
    const target = path.join(out, view + ".svg");
    await download.saveAs(target);
    const svg = fs.readFileSync(target, "utf8");
    for (const expected of [
      "10.1093/bioinformatics/btab590",
      "10.1093/bioinformatics/btw039",
      "UniProt 2026.02",
      "Distinct accepted peptides",
      "Unipept Visualizations",
    ])
      if (!svg.includes(expected))
        throw new Error("Missing export citation/context " + expected);
  }
}
await page.locator("[data-view=sunburst]").click();
await page.selectOption("#metric", "signal");
const promise = page.waitForEvent("download");
await page.locator("#save-png").click();
await (await promise).saveAs(path.join(out, "sunburst_signal_export.png"));
await page.screenshot({
  path: path.join(out, "full_report.png"),
  fullPage: true,
});
// Exercise empty data and long identifiers without changing the original report.
const document = fs.readFileSync(path.resolve(process.argv[2]), "utf8");
const match = document.match(/const data = (.*);/);
if (!match) throw new Error("Missing report payload");
const fixture = JSON.parse(match[1]);
const empty = structuredClone(fixture.samples[0]);
empty.sample = "Synthetic empty acquisition";
for (const measure of ["peptides", "signal"]) {
  empty.totals[measure] = 0;
  empty.trees[measure] = {
    id: 1,
    name: "All evidence",
    rank: "root",
    count: 0,
    selfCount: 0,
    children: [],
  };
  for (const category of Object.values(empty.balance)) category[measure] = 0;
  empty.overlap[measure] = 0;
}
fixture.samples.push(empty);
fixture.samples[0].sample = "Synthetic long label " + "W".repeat(220);
const serialized = JSON.stringify(fixture)
  .replaceAll("<", "\\u003c")
  .replaceAll("&", "\\u0026");
const edgepath = path.join(out, "edge_fixture.html");
fs.writeFileSync(
  edgepath,
  document.replace(match[0], "const data = " + serialized + ";"),
);
await page.goto("file://" + edgepath);
await page.waitForSelector("body[data-ready=true]");
await page.selectOption("#sample", String(fixture.samples.length - 1));
if (
  !(await page.locator("#save-svg").isDisabled()) ||
  !(await page.locator("#save-png").isDisabled())
)
  throw new Error("Empty evidence still offers figure export");
if (
  !(await page
    .locator("#taxonomy-stage")
    .textContent()
    .then((t) => t.includes("No measured evidence")))
)
  throw new Error("Missing empty-data explanation");
await page.selectOption("#sample", "0");
const longDownload = page.waitForEvent("download");
await page.locator("#save-svg").click();
const longFile = path.join(out, "long_label.svg");
await (await longDownload).saveAs(longFile);
const wrapped = fs.readFileSync(longFile, "utf8");
if (wrapped.includes("W".repeat(220)))
  throw new Error("Long identifier was not wrapped");
if (!wrapped.includes("10.1093/bioinformatics/btw039"))
  throw new Error("Citation lost in wrapped export");
if (errors.length || remote.length)
  throw new Error(JSON.stringify({ errors, remote }));
fs.writeFileSync(
  path.join(out, "COMPLETE.json"),
  JSON.stringify(
    {
      status: "PASS",
      views: 3,
      metrics: 2,
      samples: 2,
      errors,
      remote,
      export_citations_checked: true,
      empty_samples_checked: true,
      long_labels_checked: true,
    },
    null,
    2,
  ) + "\n",
);
await browser.close();
