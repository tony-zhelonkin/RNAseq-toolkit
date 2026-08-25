# bulkiRNA — handoff and plan of record

**Updated:** 2026-08-21 · **Branch:** `feat/bulkirna-package` · **Last release:** `v1.0.0` = `0ceda46`
**Development version:** `1.1.0.9000`, naming the *next* release rather than this tree.
**Gates, measured at the release commit:** 1,784 tests, 0 failures, 6 environmental skips ·
golden 20/20 · 7 identity assertions ·
`devtools::check(error_on = "note")` → 0 errors, 0 notes, 1 warning (`qpdf` absent from the image).
**Quote the invocation with the verdict:** `check()` returns rather than failing, so a run printing
warnings still reads as a pass without `error_on=`.
**Surface:** 59 exports — 47 stable, 12 experimental, **0 deprecated** · 3 signature-frozen · 6 stochastic
**Retained, not exported:** the 21 legacy names, as golden-baseline fixtures.

> **Do not copy gate numbers forward.** They were carried between sessions once and were stale: the
> suite was red on two tests while being reported green. Re-run before quoting.

This is the entry point. Everything else is reachable from here.

| Document | What it owns |
|---|---|
| [2026-08-10-bulkirna-package/00_INDEX.md](2026-08-10-bulkirna-package/00_INDEX.md) | Phases 0–7 and their execution record |
| [2026-08-13-analysis-api-roadmap/00_ROADMAP.md](2026-08-13-analysis-api-roadmap/00_ROADMAP.md) | Phases 8–9. **§7–§16 are the execution record**, including every gate that was superseded and why |
| [2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md](2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md) | **The canonical previous implementation of each analysis, with paths, parameters and the reasoning.** Read this instead of asking where prior work lives |
| [2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md](2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md) | What CoReSh's score is, traced to its implementation |
| [2026-08-13-analysis-api-roadmap/03_DEFERRED.md](2026-08-13-analysis-api-roadmap/03_DEFERRED.md) | TF activity, PROGENy, WGCNA — parked, with findings intact |
| [2026-08-13-analysis-api-roadmap/04_NAME_AUDIT.md](2026-08-13-analysis-api-roadmap/04_NAME_AUDIT.md) | Why each public name and argument is spelled as it is |
| [2026-08-13-analysis-api-roadmap/05_RETURN_AUDIT.md](2026-08-13-analysis-api-roadmap/05_RETURN_AUDIT.md) | What every export returns, with the reasoned exceptions |
| [2026-08-13-analysis-api-roadmap/06_CONSUMER_MIGRATION.md](2026-08-13-analysis-api-roadmap/06_CONSUMER_MIGRATION.md) | **Phase 4's recipe:** the call mapping, the four ways a migrated script fails while looking fine, and how to verify one |
| [2026-08-10-bulkirna-package/07_SUBMODULE_RETIREMENT.md](2026-08-10-bulkirna-package/07_SUBMODULE_RETIREMENT.md) | **Phase 6's procedure**, per repository, with rollback and stop points. Its measured input is `07_SUBMODULE_STATE.txt`. Not executed |
| [2026-08-11-coresh-extraction/00_PLAN.md](2026-08-11-coresh-extraction/00_PLAN.md) | CoReSh extraction, C0–C6 |
| [../adr/](../adr/) | The four architectural decisions |
| [../../../AGENTS.md](../../../AGENTS.md) | The rules an agent must follow in this repo |

---

## 0. Why this project exists

`RNAseq-toolkit` was a folder of scripts consumed as a git submodule and pulled in with
`source()`. 24 copies, each pinned to a different commit, so a fix in one reached nobody. Nothing
was tested. Nothing carried a version, so no result could name the code that produced it. **The
goal is not tidiness. It is that a number in a figure can be traced to a version of the code that
produced it, and that a defect fixed once is fixed everywhere.**

The second goal, stated by the owner, is why Phases 8–11 exist at all:

> *a thin clean API that just does the job the way I like it, without having to re-state to agents
> implementing analysis for me to go search for my previous experience reference repos — a
> zero-token architecture for my personal preferences on how bulk and pseudo-bulk should be done.*

That reframes the package as **the executable record of a set of methodological preferences**, and
it sets the acceptance bar in three ways:

1. **A function needing a paragraph of instructions to use correctly has failed.** The defaults
   *are* the preferences, documented once.
2. **A choice made the same way three times is a default, not a decision.** Three projects agreeing
   on `top = 500` for PROGENy is no longer a judgement call.
3. **Prose is not enough.** `bulkirna_api()` and `bulkirna_stochastic()` are enforced in code with
   tests; the pipeline conventions in `01_REFERENCE_PROJECTS.md` are prose and can drift the way the
   originals did. Where a preference matters it should become a default with a test, and the
   remaining prose should say plainly that it is prose.

### For an agent asked to run an analysis "the way it is usually done here"

Read in this order, and do not re-derive any of it:

1. `AGENTS.md` — the rules, the container invocation, the gates.
2. `vignettes/` — the order of operations per layer: `gene-sets`, `differential-expression`,
   `gatom`, `coresh`.
3. `01_REFERENCE_PROJECTS.md` §0 and §0b — the conventions that hold across every stage, and the
   house `1.x` ingest / `2.x` compute / `3.x` render pipeline shape with its config keys and output
   layout.
4. `01_REFERENCE_PROJECTS.md` §1–§4 — the canonical prior implementation of the specific analysis,
   with the project path, the parameters used, and where those projects disagreed with each other.
5. `bulkirna_api()` for what is safe to build on; `bulkirna_stochastic()` for what is random.

---

## 1. The pain points

The first four are original; the rest were found during the work, each measured rather than inferred.

**1. Silent wrongness that a green run cannot distinguish from a null result.** Ten confirmed:

| # | Defect | How it hid |
|---|---|---|
| 1 | `msigdbr` keys its ortholog cache on species alone | Second collection in a session truncated to the intersection; filed as issue 62 |
| 2 | A ten-column allowlist NA-filled four columns | 492 rows with `padj = 2.6e-33` beside `neg_log_padj = NA` |
| 3 | `annotate_genes(use_biomart = TRUE)` returned all-NA | Silently, after being explicitly asked for biomaRt |
| 4 | `safe_install()` logged two rows per failure | One falsely "returned without error" |
| 5 | `load_or_compute()` served results older than their input | Table listed 9 sets that no longer existed, omitted 27 that did |
| 6 | `gatom::scoreGraph()` is stochastic and was never seeded | Same input scored 165.85 then 169.32 in one session |
| 7 | A cached `NULL` made a failure outlive its own fix | Next run reads the file, skips the contrast, reports the same clean nothing |
| 8 | A green gate that never ran the code | Exit 0 printing `cache hit`; forcing a recompute showed every cell had moved |
| 9 | `bplapply()` switches `RNGkind()`, so seeding without pinning changes answers | `gs_test()` returned different p-values inside a parallel worker, no warning |
| 10 | `NA` Entrez ids in 0.7% of real CoReSh datasets | Passed the entire test suite; found by a compendium sweep |

Numbers 5–10 share one shape: **an empty result and a broken run were indistinguishable from inside
the code, so the broken run took the empty result's path.** In one case that purged 492 rows from two
master tables and exited 0.

**And the shape recurred inside the package, after being named.** `coresh_sets()` counted extractions
rather than sets, so a run where every set fell outside the size bounds returned a clean empty `gs_db`
with no message — the reference implementation's `skip:` hole moved one step downstream rather than
closed. **This is the failure mode to assume is present.**

**2. Machinery retyped per project, then drifted.** ~81% of the surveyed CoReSh lines were generic;
the engine was md5-identical in four checkouts, re-vendored by a `cp` in a comment.

**3. The same drift reproduced *inside* the package.** Seeding was implemented three times, three
ways. Then five species handlers disagreed about what a species is —
`gene_to_entrez(species = "Homo sapiens")` worked while `coresh_search(species = "Homo sapiens")`
errored, same session, same species. TF activity is next in line: two estimators and three network
provenances across three projects.

**4. Reference data fetched live, so a result depends on the day.** `msigdbr`'s cache, the GATOM
downloads, and a documented OmniPath breakage.

**5. Provenance that can be silently meaningless.** `.ref_path()` recorded
`basename(Sys.readlink(current))`; a `current` pointing outside its source, or being a real
directory, recorded a snapshot identifying nothing while the call succeeded. **And the opposite
failure:** `gsdb_coresh()`'s first version `dput()`d the queries and the whole hits table into the
provenance to satisfy a scalar validator — kilobytes of escaped R code that looked like completeness
and was unreadable. An honest gap a reader can see beats a blob that looks complete.

**6. The agent-facing documentation described the dead predecessor.** `AGENTS.md` — the one file
every assistant reads, via `CLAUDE.md` and `GEMINI.md` — was 225 lines about a script library,
including "there is no `library(RNAseqToolkit)`". Rewritten 2026-08-14.

**7. A document can record a method that provably has a gap.** `04_NAME_AUDIT.md` said its table was
"measured from the formals of all 57 live exports" — the technique that cannot see an argument
declared on an S3 method behind a generic's `...`, which is how `group` and `samples` went unaudited.
The counts were right; the counting was not. **Those are separate claims and only the first was
true.**

**8. A false premise attached to an enforcement-test exception does active work.** Two exceptions in
the argument audit claimed their spellings were "fixed by frozen signatures". None of the seventeen
exports involved was frozen and none had a consumer call site. A wrong row in a plan document is
merely wrong; **a wrong reason in a test holds open the gap the test exists to close**, under a
sentence that reads as settled.

**9. An enforcement test's success path is untested until the day it succeeds.** Removing the last
multi-spelling exception turned the suite red: the comparison had only ever run with a non-empty
allowlist, and at zero the two sides had different shapes. The same inversion is scheduled for
`v1.0.0`, when the 21 shims go and every shim allowlist becomes permanently empty. Checked early and
fixed; the helpers now carry "do not delete as dead code" and the reason.

**10. A stale installed copy of the package makes every downstream claim unverifiable.** Phase 4
resumed against `bulkiRNA 0.3.1` in the shared library, which predates `gs_to_master()` — the exact
verb the migration targets. The development tree said 0.5.0.9000 and all the gates were green, so
nothing in the repo pointed at the problem. **The library a consumer loads is a fourth version
authority, and ADR-001 did not cover it.** Now the first precondition in
[06_CONSUMER_MIGRATION.md](2026-08-13-analysis-api-roadmap/06_CONSUMER_MIGRATION.md) §1.

**11. A registry can name a successor that does not cover the call site.** `bulkirna_api()` gives
`convert_human_to_mouse`'s successor as `gsdb_msigdb(species = ..., db_species = "HS")`, which is
right for MSigDB and only for MSigDB. The real consumer calls it on a T2G parsed from a custom GMX,
and no export ortholog-maps an arbitrary `gs_db`. The registry is machine-readable, tested, and
non-empty, and was still wrong for this use — **the deprecation metadata is tested for presence and
shape, not for whether the successor can do the job.** Found by migrating, which is the only way it
could have been found.

**12. The shape recurred a third time, in a consumer's own error handling.** The script being migrated
wrapped its normalisation in `tryCatch(..., error = function(e) NULL)`. With the §3.1 defect above,
that turned a hard type error into zero rows and no message — the migration's first run reported
`Normalized rows: 0` and carried on. `gs_validate_master()` is what stopped it, several steps later.

**13. Two enforcement tests had never executed.** `test-error-audit.R` walked the parse tree and passed
the empty symbol into its own recursion; a formal declared without a default *is* that symbol, so it
aborted on the first function it met. Both it and a newly written test called `expect_gt(info = )`,
which testthat rejects as an unused argument. Neither had ever run, and the suite reported them as
errors rather than failures — which reads past easily in a summary line. **The error audit passed on
its first real run, so what it checks was already sound; only the check was broken.** Pain point #9's
shape, live in two places while the gates were being described as green.

**14. A version string that identifies two different code states.** `DESCRIPTION` read `0.5.0` while
`HEAD` stood 50 commits and ~3,032 changed lines of `R/` past the `v0.5.0` tag, adding nine files and
15 exports. ADR-001 forbade a version moving without a tag and never forbade code moving under a
released version. Closed by [ADR-001's amendment](../adr/ADR-001-reproducibility-unit.md) and
`test-version-identity.R`. **The general asymmetry is the durable lesson: every *content* invariant
here fails a test when violated, while identity and deployment invariants were prose, because they
live between repositories, tags, builds and installed libraries where package tests do not reach.**

**15. Declared is not verified, on both sides of the image.** Seven R packages present in v0.5.10 and
absent from v0.5.13 were never declared in the build at all — transitive arrivals that left when an
upstream `DESCRIPTION` changed, with no commit recording the loss. `torch` was never declared either,
so the image's GPU capability was whatever PyPI shipped on build day: both images landed
`2.13.0+cu130`, needing driver ≥ 580 against this host's 565.57, and CUDA initialisation failed while
`nvidia-smi` still showed the GPU. **`pip` and `install.packages()` both exit 0 after resolving a
version nobody asked for**, so a pin is a request until something checks the result.

**16. A briefing artifact can be silently truncated.** The signature reference handed to a delegated
agent was generated with `deparse(args(f))[1]` — first line only — so every multi-line signature was
cut mid-formal. The agent stopped and said so rather than guessing, which is the behaviour to want.
The same artifact also propagated a wrong claim from a plan document: that `gsdb_from_file()` takes
`database_label` for the join key, when `database` is the key and `database_label` is only a renderer
display string.

**17. A process-matching pattern can capture someone else's work.** `pgrep -f "codex exec"` matched a
concurrent session's agent in a different project, so a delegated run that had died instantly was
reported as running, twice. A `pkill` on the same pattern could have killed that session's work.
**Match on the working directory, not the command name.**

**18. An omitted argument is not an absent default.** The toolkit's `run_gsea()` passed no size
bounds, and its header said "No minSize/maxSize exposed" — which reads as *no filtering* and meant
*`clusterProfiler`'s 10 and 500*. `gs_test()` has no such default, so a faithful argument-for-argument
translation would have widened every gene-set universe and shifted every padj through the
multiple-testing correction, with nothing in the diff to show it. **Check the callee's defaults for
every argument the old call omitted**, not just the arguments it passed.

**19. A count is not a membership check.** Three project data files were replaced by bundled
providers, and the migrated script verified each by asserting the provider's set count. TransportDB
gave 49 against 49 with **zero shared names**: a `.gmx` is column-oriented and had been read row-wise,
the provider keys sets by dot-joined hierarchy path where the file names only the leaf, and one
comparison was human symbols against mouse. All three errors are invisible to a count. Corrected, all
three substitutions are sound — but the assertion that was in the code would have passed either way.
**And the baseline is what the legacy path produced, not the file it read**: the MitoPathways file is
human and the old code converted it, so the raw file was never what the analysis saw.

**20. A file you must edit is not a file you may commit.** Both DC-nexus `docker-compose.yml` files
needed the image moved forward for the migrated scripts to run, and both already held another actor's
uncommitted work — their own image bump, changed ports, changed data mounts. The edits were made and
left uncommitted, because committing them would have swept someone else's in-progress changes into my
commit. `14616-DM`'s equivalent file was clean, so there the bump is committed. Same action, opposite
treatment, decided by what else was in the file. **The same rule reappeared one level down at the
submodule cleanup:** six repositories have parent gitlinks that already differ from their recorded
value, so advancing one would absorb an unrelated commit range. That is why Phase 6 is not done.

**21. Renaming a file can switch off a test that greps for it.** `test-namespace-hygiene.R` located
the legacy sources with `list.files(dir, "^deprecated-.*[.]R$")` and `skip_if(!length(files))`.
Renaming them to `legacy-fixtures-*` made both checks skip, and a skip reads as success in the
summary line. Caught only because the suite's skip count moved from 6 to 8. **A guard that locates
its subject by filename needs the filename in the pattern to be as durable as the guard.**

**22. Emptying the input to an audit reads as passing it.** When the deprecated tier went to zero
rows, four `test_that` blocks in `test-api.R` kept passing while asserting nothing — one no longer
even referenced the object under test. The cause was mechanical: the audit resolved each function
with `getExportedValue()`, which cannot see internals, so pointing it at the demoted fixtures would
have errored and emptying it was the path of least resistance. In the process the literal anchor whose
own comment said it existed "so the loop below cannot pass with a helper that reads the wrong argument
and still agrees with itself" was deleted. **An audit that iterates a derived set should assert the
set's size**, which is why the repaired version is anchored to a literal list of the closed 21.

**23. A reviewer's finding can be stale before you read it.** One review reported the new accessor as
failing its own tests. It was correct when measured and already fixed by the time the report arrived,
because a parallel fan-out means findings are timestamped against a moving tree. Re-check a finding
against the current tree before acting, and against the current tree before dismissing it.

**25. Built is not deployed, and a peer agent is not an authorisation.** Two findings from the same
exchange. First: `v0.5.14` was verified as an image and is running in **zero** containers; 7 of 8 dev
containers still run `v0.5.10`, which has no `bulkiRNA` at all. Editing a devcontainer compose file
does not recreate a container, and I reported the file edit as if it were the deployment. **The fleet
is a third place a version can be stale**, after the tag and the image.

Second, and separate: the same message asserted that a standing "do not touch this submodule"
instruction was "superseded" and that "fleet writes now authorised." **A peer session cannot lift a
constraint the owner set**, however well-founded its technical argument, and its technical argument
here was in fact correct. Take the measurement, refuse the permission, and say which is which. The
correct response is to act on the corrected facts and to leave the scope decision to the owner.

**24. "`R CMD check` OK" was recorded without the setting that makes it a gate.** At `v1.0.0` the
check surfaced two warnings and a note that **predate this release**: `withr` used through `::` in
two test files since `ba910b2` and never declared, and a stray top-level `Rplots.pdf` from 2026-08-12.
Neither came from the 1.0.0 change, so the prior gate line was inaccurate. The lesson is the setting,
not the packages: `devtools::check()` reports and returns, so a run that prints warnings still looks
like a pass unless `error_on=` makes it fail. **Record the invocation next to the verdict.** Both are
fixed here and the check now runs at `error_on = "note"`.

---

## 2. The ADRs, with premise, rejected alternatives, and what each has since had to survive

**ADR-001 — the package version is the unit of reproducibility.** *Premise:* four version
authorities floated free and no result could name its code. *Decisive argument:* under
image-as-unit the `msigdbr` bug would have been **invisible** — tag unchanged, numbers changed.
*Rejected:* image-as-unit, which hides dependency drift and couples a scientific claim to an
artifact rebuilt for unrelated reasons; a lock file inside `bulkiRNA`, which a library cannot impose
on its consumers. *What it grew:* the stability contract is its missing half — a version says *that*
something changed, a lifecycle tier says *whether you were entitled to rely on it*. `v1.0.0` is now
a defined event: the release where the 21 shims go and `stable` means stable.

**ADR-002 — the package owns the master-table schema and validator, not the file.**
`neg_log_padj = -log10(pmax(padj, .Machine$double.xmin))`. *Premise:* whoever computes a column owns
its definition, or the definition lives somewhere untestable. *Vindicated:* `37.55` in the live
table, above the retired cap of 16, on real values. *Rejected:* toolkit ownership — a skill cannot be
`library()`d, so the schema reverts to prose, which is how it drifted three ways; per-project schema
files, which drift by default. *Held under pressure twice:* when `gs_test()` gained `log2err`,
`gs_to_master()` still returned exactly 14 columns; and when `gs_coregulation()` introduced an
unsigned statistic, the NES guard was left **untouched** — a `pct_var` result reaches the master table
only with `stat_as_nes = TRUE`, and a test pins that the guard fires by default. A master table whose
`nes` column silently holds three different statistics is the defect this project began with.

**ADR-003 — `Suggests` stays.** *Premise:* it is the machine-readable list `R CMD check`,
`install_deps` and the image build read. *Decisive argument:* 14 light `Imports` are why
`install_github` works on bare R in seconds. *Rejected:* collapsing into `Imports`, making every
optional feature a hard install blocker; core-plus-feature packages, cleanest isolation but a second
release train for one maintainer. *Amended twice on `Remotes:`* — unnecessary, then necessary when
`coresh` looked GitHub-only, then unnecessary again when `coresh` turned out to export nothing.
*Its seam proved itself:* adding `knitr` and `rmarkdown` for S6 immediately failed the
registry-versus-`Suggests` drift test, which is what that test is for. They are recorded as
**dev-only beside `testthat`**, because `bulkirna_check_deps()` answers "what must I install to use
this feature" and that answer never includes knitr. *Watch item:* WGCNA's dependency tree is the
first that makes the core/feature split arguable.

**ADR-004 — two reference-data tiers, bundled or refcache.** *Premise:* three mechanisms had grown
and a fourth environment variable was about to appear. Bundling MitoCarta stands: a source with no
API cannot be re-fetched reproducibly, so vendoring **is** the archive, and 168 KB makes the size
objection moot. *The insight:* nobody needs to own pinning, someone needs to own **recording**.
*Earning its keep three times:* CoReSh chunks resolve through it; the OmniPath networks are exactly
its case, fetchable in principle and unfetchable in practice at current versions; and `gsdb_coresh()`
now carries `snapshot = syn66227307_20260721` in the returned object's provenance rather than in a
log line, which is what the ADR asked for and had not yet been made to do.

**ADR-005 was proposed and withdrawn.** It was to settle ULM versus MLM for TF activity. With that
layer parked there is no decision to record, and an ADR for unscheduled work rots. It survives as an
open question in `03_DEFERRED.md`.

---

## 3. Where we are

| Phase | Status |
|---|---|
| 0 Inventory & freeze | ✅ |
| 1 Skeleton, internal `source()` removed | ✅ |
| 2 Prove the dev loop | ✅ against real data |
| 3 Install into the image | 🟡 `v0.5.14` verified as an image (`0.6.0` at `e42c2de`) but **running nowhere**; 7 of 8 dev containers are on `v0.5.10`, which has no bulkiRNA. Built ≠ deployed |
| **4 Migrate heavy consumers** | ✅ **complete** — `14839-DM-cGAS`, Meta-Aging (`d6633dd`), `DC_Dictionary` (`86d78cf`), `DC_hum_verse` (`d7533c0`) · STING-JR excluded by the owner. Scope and recipe in [06_CONSUMER_MIGRATION.md](2026-08-13-analysis-api-roadmap/06_CONSUMER_MIGRATION.md) |
| 4c CoReSh extraction | ✅ C0–C4 (C4 landed as G3) · C5/C6 🔒 **owner-gated**, not blocked on work |
| 5 Bind the skills | 🔒 **owner-gated**, not blocked on work · largely the same work as C6 |
| 6 Retire the submodule | 🟡 procedure written and reviewed · blocked on six dirty gitlinks, not on work |
| 7 Distribution | ⬜ · `v1.0.0` is tagged, so the surface it would distribute is now fixed |
| **8 Stabilize the surface** | ✅ **S1–S6 complete** |
| **9 One CoReSh/GESECA run** | ✅ **G1–G4** · 🚫 **G5's web-UI half needs a browser** |
| 10–11 Activity layer, WGCNA, uniform surface | 🚫 parked by decision, findings intact |

### What Phases 8–9 delivered

- **`bulkirna_api()`** — 79 exports with lifecycle and signature-freeze as *independent* axes,
  because 20 of the 24 frozen names are also deprecated and one column would have to lie about one or
  the other. Its load-bearing test compares the registry against `NAMESPACE` in both directions.
- **`bulkirna_stochastic()`** — the 7 functions that consume randomness, each seed argument, default
  and source. **The three default seeds deliberately disagree** — `coresh_*` uses upstream's literal
  `1L`, `gatom_module()` uses `42`, the fgsea adapter uses `123L` — because each matches numbers
  already published. Documented deviations, not untidiness.
- **`R/rng.R`** — the only place allowed to mutate RNG state, with a parse-tree test that fails if a
  fourth site appears, verified in both directions.
- **One species resolver** where five disagreed, widening every alias and narrowing none except
  `annotate_genes(species = NULL)`, which used to mean mouse silently.
- **One spelling per concept**: `top`/`top_n`/`n_top` → `top_n`, `color_palette`/`colours`/`palette`
  → `palette`, decided by the owner once the cost was measured at zero.
- **CoReSh, complete** — `coresh_chunks/match/search/convergence/validate/loadings/sets` plus
  `gsdb_coresh()`. `pct_var` bit-identical to hand computation on 12 of 12 real datasets; the loading
  projection identical to the reference implementation on **63 of 63** comparable hits; the provider
  identical to its own hand composition, provenance included. Dataset identity is the `(gse, gpl)`
  pair: **1,635 of 42,465** human accessions appear on more than one platform, and on the one
  measured the two platforms shared **19 of 50** top-loading genes.
- **`gs_coregulation()`** — GESECA on any expression matrix, `stat_type = "pct_var"`, `direction`
  `NA` because the statistic is unsigned.
- **Four vignettes**, two genuinely evaluated, and `R CMD check` at 0/0/0 with vignettes rebuilding.
- **A real-data fixture** carrying the three structures that have broken this package, with a
  generating script that reproduces it byte-for-byte.

### The strongest evidence the port is correct

`HALLMARK_HYPOXIA` against the whole human compendium: **44,253 datasets in 87 seconds**, and **7 of
the top 10 are explicitly hypoxia or HIF experiments**, with VHL — the canonical HIF degradation
pathway — at rank 4. Validated against biology, not against itself.

### Gates not met as written, recorded as superseded rather than quietly reinterpreted

| Gate | Why not met | What replaced it |
|---|---|---|
| G2: memberships match DC-nexus's stored GMT | 28 of 58 match; the chunks that produced that GMT no longer exist, established four ways | Loading projection identical to the reference on 63 of 63 hits |
| G4: a golden baseline added | That harness exists to prove **frozen legacy names** did not drift; a new function would also pin fgsea's estimator version into it | Agreement with the hand-computed GESECA formula, plus planted-signal and seed-relation tests |
| S6: vignettes runnable on the shipped fixture | `tests/fixtures/` and `data/` are both `.Rbuildignore`d | Gene sets on real `inst/extdata` sets; DE on declared synthetic data; GATOM and CoReSh honest and unevaluated |
| §5: run G5 before G2–G4 | G5's remaining half needs a browser | Gate each step against the reference implementation instead; G5 stays the only external falsifier |

---

## 4. Next immediate steps

**Both jobs that were in flight have landed and been verified.** `scdock-r-dev:v0.5.14` is built and
checked: `torch.cuda.is_available()` is `TRUE` under `--gpus all` with `2.11.0+cu128`, `bulkiRNA` is
`0.6.0` at `RemoteSha e42c2de`, all nine required packages load, 79 exports, and `verify_base.py`
passed inside the build across 24 pins. `scrublet` and `skimage` import from `/opt/venvs/base`.
`14616-DM` is on `v0.5.14` (commit `16e79bc`).

**Phase 4 is complete.** All four consumer groups are migrated; only the owner's STING-JR exclusion
remains outstanding as a decision, not as work.

1. **The `v1.0.0` shim decision is now the critical path.** Phase 4 is done, so nothing further will
   reduce the 21 shims: they exist because STING-JR still calls the old names, and it is excluded by
   the owner's decision (survey kept in §0 of
   [06_CONSUMER_MIGRATION.md](2026-08-13-analysis-api-roadmap/06_CONSUMER_MIGRATION.md)). The choice
   is to migrate STING-JR after all, or ship `v1.0.0` with the shims present — and the second
   **reopens S4**, which recorded all 21 as removed in `1.0.0`.
2. **Close the two gaps Phase 4 exposed**, both small and both `v1.0.0` blockers:
   an ortholog verb for arbitrary `gs_db` objects, or a narrowed `superseded_by` for
   `convert_human_to_mouse`; and an exported accessor for the master-table columns, so consumers stop
   reading `inst/extdata/master-schema-v1.csv` directly — whose **row order is not the column order**
   `gs_validate_master()` requires, which is a trap in itself.
3. **The RSPM snapshot pin does not apply to the CRAN stack.** `install_core.R:7` sets the
   `2026-04-15` snapshot as the default `repos`, and then **all 11 install calls that matter override
   it** with `repos = "https://cloud.r-project.org"` (9 in `install_core.R`, 2 in `install_httpgd.R`).
   So the pin that looks like the determinism mechanism is bypassed. This sits beside the other dead
   mechanism in the same file: `install_renv_project.R` restores `/opt/settings/renv.lock` when
   present, and the Dockerfile copies no lockfile there, so that branch never runs. `docs/build.md`
   no longer claims determinism (`scbio-docker` `3dd9b49`), but the installer is unchanged.
   Bioconductor is constrained to release `3.22`, not to versions; three r-universe sources float;
   GitHub installs float except `bulkiRNA`, which is commit-pinned. **Removing the overrides is the
   cheap fix and needs a build to verify** — they may exist because the snapshot lacks some package.
   **This is the mechanism behind pain point #15.**
4. **G5's other half** — the same query through <https://alserglab.wustl.edu/coresh>, compared
   accession by accession. Needs a browser, so it is the owner's step, and it is the only remaining
   independent check on the CoReSh port.
5. **Phase 6 — blocked on six dirty gitlinks, not on work.** Zero `source()` calls and zero legacy
   function calls reach the vendored toolkit in either project; nine checkouts remain as dead weight.
   A reviewed removal procedure exists at `/tmp/wsC/PROCEDURE.md` — **copy it somewhere durable before
   `/tmp` is cleared.** It stops short of executing because six repositories already have parent
   gitlinks differing from their recorded value: `14616-DM`, `14761-DM`, `14782-DM`, `DC_Dictionary`,
   `DC_hum_verse`, `DC_mouse_cancer`. Advancing any of them absorbs an unrelated commit range —
   `14616-DM`'s pending range mixes the bulkiRNA migration with someone else's catalogue work, and
   `14782-DM` is on another actor's `feat/consensus-package` branch. **The owner must commit or stash
   in those six first.** Two path spellings exist (`01_modules/RNAseq-toolkit` and, inside
   `DC_Dictionary`, `01_scripts/RNAseq-toolkit`), all at `752481f` — the same commit the golden
   baseline was captured at. `01_modules/SciAgent-toolkit` is adjacent in every relevant
   `.gitmodules` and is out of bounds, so any edit there must be surgical.
6. **Then reconsider Phases 10–11** with `03_DEFERRED.md` in hand: TF activity, PROGENy, WGCNA. The
   sequencing argument for parking them was that a default in a package is load-bearing in a way a
   default in a script is not, and those three are where the methodological drift is worst. That
   argument is unchanged; what has changed is that the surface they would sit on has stopped moving.

---

## 4b. The three skills, and why they are owner-gated rather than blocked

Phase 5 and CoReSh C5/C6 are the same body of work seen from two plans. All of it lives in
`scbio-docker/toolkits/SciAgent-toolkit/skills/`, measured 2026-08-21:

| Item | Skill | Current state |
|---|---|---|
| C5 | `coresh-signature-search` | Still ships `scripts/{coresh_batch,extract_gene_loadings,symbols_to_entrez}.R`, sourced into one another. That logic is now `gsdb_coresh()` and the `coresh_*` verbs. |
| C6 / Phase 5 | `bulk-rnaseq-gsea` | `references/msigdb.md` describes "the RNAseq-toolkit wrapper around `clusterProfiler::GSEA()`"; `references/visualization.md` cites `01_scripts/RNAseq-toolkit/scripts/GSEA/…` paths. |
| C6 / Phase 5 | `annotate-bulk-rnaseq-data` | `SKILL.md:71` pins **RNAseq-toolkit v0.2.0** by name; `references/te-annotation.md` also pins TE-RNAseq-toolkit v0.1.0. |

**Corrected 2026-08-24, after a challenge from the SciAgent-toolkit session.** I wrote that these
were gated only on owner sequencing, quoting `00_PLAN.md` §8.5's "what remains is sequencing against
the other track, not capability." **That quote is true of an image and false of the fleet.** Measured:

| Image | bulkiRNA | Running dev containers |
|---|---|---|
| `scdock-r-dev:v0.5.10` | **absent** | **7** |
| `scdock-r-dev:v0.5.13` | 0.4.0, 64 exports | 1 (`jr-mc`) |
| `scdock-r-dev:v0.5.14` | 0.6.0 at `e42c2de` | **0** |

So **no running container can reach a current `bulkiRNA`**, and 0.4.0 exports **no `coresh_*` verb at
all** — the CoReSh layer landed after that tag. Rewriting C5 to name `coresh_search()` would replace a
stale instruction with an unrunnable one. The other twelve verbs the migrated consumers call *do* exist
in 0.4.0, so those scripts would execute there; but they were written and gated against `1.0.0`, and
producing results from 0.4.0 under code gated at 1.0.0 is precisely what ADR-001 exists to forbid.

**So there are two gates, not one.** The sequencing decision is the owner's. The delivery gap is work:
containers must be recreated onto an image carrying the current package. That second gate is real and
was mine to notice.

**And a claim of mine to retract:** I reported "`14616-DM` moved onto `v0.5.14` (`16e79bc`)". The
compose *file* says `v0.5.14`; the *running container* is still `v0.5.10`. Editing a compose file is
not recreating a container. That is pain point #15's shape — declared is not deployed — committed by
me, in the same session that named it.

**Why it is worth doing despite being small.** §0's goal is that an agent need not be re-told these
preferences, and a skill is the first thing an agent reads. Today two of these three will send it to
`clusterProfiler` and to submodule paths no consumer uses, and one will pin it to `v0.2.0` while the
package is at `v1.0.0`. The code no longer misdescribes itself; these three files still do.

---

## 5. Open items that block nothing

- The image still pins `bulkiRNA@v0.4.0` (`scbio-docker/docker/base/R/install_core.R:174`), so nobody
  using `scdock-r-dev:v0.5.13` has the GATOM reproducibility fix. A 150-minute rebuild is the owner's
  judgement call.
- `qs2` is absent from the image; every CoReSh chunk reader needs it, and since it is a declared
  `Suggests`, `R CMD check` needs a scratch library mounted at `R_LIBS` until the image carries it.
- The consumer GMT that G2 was gated against is **not reproducible from any snapshot on disk**. Its
  companion `coresh_provenance.csv` is also wrong — 13 distinct `pctVar` values across 58 rows, one
  per query rather than per hit.
- `gatom_download_refs()` still defaults to the repo-relative `00_data/references/gatom`. The frozen
  shim must keep it, but a downloader whose default destination is relative to the working directory
  is how the test suite once wrote 16 MB into `tests/testthat/`.
- **Both source-tree enforcement tests skip under `R CMD check`** — the RNG-ownership parse-tree test
  and the name audit — because an installed package has no `.R` files. A green check does not cover
  them; only `devtools::test()` does, and the skips say so.
- The argument-vocabulary test is **157 of 178 formals in two generic buckets**. It cannot rot, and it
  enforces snake_case with reasoned exceptions, but a new spelling can still be waved through as a
  singleton concept. Stated plainly rather than dressed up.
- `10_gatom_modules.R`'s combined-network path still calls `gatom::` directly.
- `14839-DM-cGAS` is uncommitted, with the moved GATOM numbers in it. Committing in a tree holding
  live research data is the owner's call.
- The `msigdbr` re-run's 47,988 → 72,408 MSigDB rows are still unread — the largest unexamined result
  in the project.
- Two `scbio-docker` commit messages say "the integrator's"; fixing means force-pushing published
  history.

---

## 6. How this work has been run, and what it cost

Implementation is delegated to `codex gpt-5.6-sol` in isolated clones; review to independent Opus
agents; **gating and merging stay with the integrator**, because the delegated agents have no Docker
socket and therefore cannot run a single gate. Every "tests pass" in this document was run here.

Four things are worth carrying forward.

**Agents stopping on a false premise is the highest-value behaviour in the loop.** Twelve times,
right twelve times: an off-by-one export count; five shims with no machine-readable successor; three
formals assumed to exist; an undecided RNG policy; an inconsistency between two briefs; a clone made
from a commit predating the work; a `gs_db` provenance premise that was false in three ways; a
`stat_type` vocabulary that would have had to be invented; a golden-case registration my own
read-only rule forbade; a fifth species handler I had missed; and vignette data that does not exist
in a built package. Briefs say "if a premise is false, stop and tell me", and should keep saying it.

**A review finding is not automatically a defect, and neither is a reviewer's method.** One review's
most severe finding argued from R's C sources that a bare `RNGkind()` initialises `.Random.seed`; one
command showed it false in R 4.5.3. Its second half — that the branch had no test — was right and
more useful. Later, checking a reviewer's own verification method found two public formals that
nobody's audit had ever seen. **Check before acting; the check is usually cheap.**

**Measure the claim, including your own prose.** Four findings across two reviews surfaced only
because a sentence was executed instead of read: the 88% singleton figure against a stated threshold,
seventeen non-frozen exports behind a "frozen signatures" exception, two unaudited S3 method formals,
and an empty-versus-empty comparison that had never run. The findings that took longest to surface
were the ones where a sentence sounded true.

**`devtools::test()` is not the gate.** It missed a semantic merge conflict between two branches that
each passed alone, stale `man/` twice, a fixture-backed test that must skip in a built package, and
an error-message rewrite that broke a matcher. **Run `R CMD check` on the built package after every
merge**, not only after an edit.
