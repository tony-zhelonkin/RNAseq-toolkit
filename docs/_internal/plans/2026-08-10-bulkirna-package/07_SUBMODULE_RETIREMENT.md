# Procedure: retire the vendored `RNAseq-toolkit` submodules

## Scope and non-negotiable rules

This is an execution procedure derived only from `reference/SUBMODULE_STATE.txt`. It has not been run against either target.

- Remove only the exact `RNAseq-toolkit` declaration named in the inventory. Never use `git submodule deinit --all`, `git submodule foreach`, `git add .`, `git commit -a`, or a recursive removal of `01_modules` or `01_scripts`.
- Do not edit, deinitialise, stage, or commit `01_modules/SciAgent-toolkit`. Record its `.gitmodules` values before each transaction and prove they are identical afterward.
- Do not use `-f` with `git submodule deinit`, `git rm`, or `git push`. A refusal is a safety stop, not an error to override.
- Do not rewrite published history. In particular, do not rebase, filter, amend a published commit, or force-push as part of this cleanup.
- Run every command below individually. Stop at each stated stop point and review the output before continuing.

## Exact inventory and order

The two roots are independent; their sequences may be scheduled independently. Within each root, finish and publish every child-repository removal before updating the parent gitlinks or removing the root's own copy.

| Order | Repository | Declaration type | Exact submodule path/name | Measured RNA branch key | Measured parent state relevant to ordering |
|---:|---|---|---|---|---|
| M1 | `/data1/users/antonz/projects/Meta-Aging/14616-DM` | sub-repository entry | `01_modules/RNAseq-toolkit` | absent | `14616-DM` is already modified in `Meta-Aging` |
| M2 | `/data1/users/antonz/projects/Meta-Aging/14761-DM` | sub-repository entry | `01_modules/RNAseq-toolkit` | `dev` | `14761-DM` is already modified in `Meta-Aging` |
| M3 | `/data1/users/antonz/projects/Meta-Aging/14782-DM` | sub-repository entry | `01_modules/RNAseq-toolkit` | `dev` | `14782-DM` is already modified in `Meta-Aging` |
| M4 | `/data1/users/antonz/projects/Meta-Aging/neuroimmune-receptor-atlas` | sub-repository entry | `01_modules/RNAseq-toolkit` | `dev` | parent gitlink is not listed as modified |
| M5 | `/data1/users/antonz/projects/Meta-Aging` | superproject entry | `01_modules/RNAseq-toolkit` | `dev` | root branch is `main` |
| D1 | `/scratch/current/antonz/projects/DC-nexus/DC_Dictionary` | sub-repository entry | `01_scripts/RNAseq-toolkit` | absent | `DC_Dictionary` is already modified in `DC-nexus` |
| D2 | `/scratch/current/antonz/projects/DC-nexus/DC_hum_verse` | sub-repository entry | `01_modules/RNAseq-toolkit` | absent | `DC_hum_verse` is already modified in `DC-nexus` |
| D3 | `/scratch/current/antonz/projects/DC-nexus/DC_mouse_cancer` | sub-repository entry | `01_modules/RNAseq-toolkit` | `dev` | `DC_mouse_cancer` is already modified in `DC-nexus` |
| D4 | `/scratch/current/antonz/projects/DC-nexus` | superproject entry | `01_modules/RNAseq-toolkit` | `dev` | root branch is `main`; unrelated tracked and untracked work is present |

An absent branch key does **not** mean `main`; it means `.gitmodules` does not select a branch. The configured submodule branch is also not evidence of the child repository's currently checked-out branch.

## Owner hold before touching six child repositories

Do not begin M1, M2, M3, D1, D2, or D3 yet. Each is already reported as a modified gitlink in its parent. Before this procedure changes that child, its owner must review these commands from the parent:

```bash
git status --short -- <child-path>
git diff --submodule=log -- <child-path>
git -C <child-path> status --short
git -C <child-path> branch --show-current
git -C <child-path> rev-parse HEAD
```

The owner must decide whether the checked-out child HEAD is the intended starting point and whether any internal work must be committed or stashed. If it is only a pointer mismatch, the owner must either approve carrying that complete commit range into the later parent gitlink update or reconcile it in a separate owner-controlled change first. There is no path-level staging technique that can split one gitlink transition into “pre-existing” and “cleanup” portions.

M4 has no measured parent pointer mismatch, but its internal status and current branch were not included in the evidence. It must still pass the gates below. The same is true of all seven child repositories.

## Phase 1: dry-run and safety gate in each repository

Use the following transaction once for every inventory row, with the literal values assigned in the invocation list later in this document. For child repositories, the owner hold above comes first.

Set and inspect the exact target:

```bash
repo_root='<literal repository path from the invocation list>'
sub_path='<literal submodule path from the invocation list>'
sub_name="$sub_path"
cd "$repo_root"
git rev-parse --show-toplevel
git status --short
git branch --show-current
git ls-files --stage -- "$sub_path"
git config -f .gitmodules --get "submodule.${sub_name}.path"
git config -f .gitmodules --get "submodule.${sub_name}.url"
git config -f .gitmodules --get "submodule.${sub_name}.branch"
```

Effects: the assignments select one literal inventory row; `cd` enters it; the remaining commands display the actual repository root, dirty state, current branch, `160000` gitlink, URL, and optional branch without changing them.

Stop unless all of these are true:

- `git rev-parse --show-toplevel` is exactly `$repo_root`.
- The current repository branch is a named, owner-approved branch. For M5 and D4 it must still be `main`; do not infer a child branch from `.gitmodules`.
- `git ls-files --stage` returns exactly one entry with mode `160000`.
- The configured path and URL match the inventory and `git@github.com:tony-zhelonkin/RNAseq-toolkit.git`.
- The branch command prints `dev` only for the rows marked `dev`; exit status 1 with no output is expected for rows marked absent.

Prove that the index and the two paths this transaction will touch are not already carrying somebody else's change:

```bash
git diff --cached --quiet
git diff --quiet -- .gitmodules
git diff --cached --quiet -- .gitmodules
git diff --quiet -- "$sub_path"
git diff --cached --quiet -- "$sub_path"
```

Effect: each command is read-only and exits zero only when the relevant comparison is clean. Stop if any exits nonzero. An empty index is required so the later commit cannot absorb pre-existing staged work. Unrelated **unstaged** files may remain, provided neither `.gitmodules` nor the exact target path overlaps them.

Check the checkout itself and prove its commit is recorded and recoverable:

```bash
gitlink_sha=$(git ls-files --stage -- "$sub_path" | awk '$1 == "160000" {print $2}')
checkout_sha=$(git -C "$sub_path" rev-parse HEAD)
test "$gitlink_sha" = "$checkout_sha"
git -C "$sub_path" status --porcelain=v1
git -C "$sub_path" remote get-url origin
git -C "$sub_path" fetch --tags origin
remote_refs=$(git -C "$sub_path" for-each-ref --format='%(refname)' --contains "$checkout_sha" refs/remotes/origin/ refs/tags/)
test -n "$remote_refs"
printf '%s\n' "$remote_refs"
git -C "$sub_path" for-each-ref --format='%(refname) %(objectname) %(upstream:short) %(upstream:track)' refs/heads/ refs/tags/
git -C "$sub_path" fsck --full --no-reflogs --unreachable
```

Effects: the first three commands prove the checkout is exactly the commit recorded by the parent; `status` must print nothing; the remote commands identify and refresh `origin`; and the containment commands prove at least one fetched remote-tracking branch or tag contains the checked-out commit. The final two read-only commands inventory local branch/tag tips and unreachable objects. Stop on a mismatch, any status output, a missing/unapproved remote, a failed fetch, an empty remote-ref result, or any local ref/object the owner has not either preserved or explicitly judged disposable.

Search every file outside the checkout, including untracked files:

```bash
rg -n --hidden --glob '!.git/**' --glob '!.gitmodules' --glob "!${sub_path}/**" 'RNAseq-toolkit|01_(modules|scripts)/RNAseq-toolkit' .
rg -n --hidden --glob '!.git/**' --glob '!.gitmodules' --glob "!${sub_path}/**" 'source[[:space:]]*\(' .
```

Effect: the first command finds any explicit outside reference to the toolkit or its paths; no output (ripgrep exit 1) is the required result. The second lists every remaining `source()` call outside the checkout; review every result and stop if any target the toolkit directly or through a constructed path. Also stop if project-specific checks reveal any live use of former toolkit functions.

Record the exact out-of-bounds stanza and local config before mutation:

```bash
sciagent_before=$(git config -f .gitmodules --get-regexp '^submodule\.01_modules/SciAgent-toolkit\.')
toolkit_local_before=$(git config --local --get-regexp "^submodule\.${sub_name}\." || true)
printf '%s\n' "$sciagent_before"
printf '%s\n' "$toolkit_local_before"
```

Effect: these commands capture the SciAgent values for the later equality check and record the toolkit's local config keys for rollback. Preserve the displayed output in the reviewed work record.

**STOP POINT A:** do not run Phase 2 until the owner has reviewed all Phase 1 output, all tests pass, both searches are safe, the target checkout is clean and remotely recoverable, and the current branch is approved.

## Phase 2: remove one declaration and quarantine its administrative repository

Resolve and validate the administrative git directory before changing anything:

```bash
module_gitdir=$(git rev-parse --git-path "modules/$sub_path")
retired_gitdir="${module_gitdir}.retired-rnaseq"
printf '%s\n' "$module_gitdir" "$retired_gitdir"
test -d "$module_gitdir"
test ! -e "$retired_gitdir"
```

Effect: these commands resolve the correct `.git/modules/...` location even when the current repository is itself a submodule, and ensure the quarantine destination is unused.

Run the removal in this order:

```bash
git submodule deinit -- "$sub_path"
git config -f .gitmodules --remove-section "submodule.${sub_name}"
git add -- .gitmodules
git rm --cached -- "$sub_path"
if test -d "$sub_path"; then rmdir -- "$sub_path"; fi
if git config --local --get-regexp "^submodule\.${sub_name}\."; then git config --local --remove-section "submodule.${sub_name}"; fi
mv -- "$module_gitdir" "$retired_gitdir"
```

Effects, in order:

1. `deinit` removes the clean working-tree checkout and normally removes its exact local `.git/config` section; without `-f`, it stops rather than discarding detected work.
2. `git config -f` surgically removes only the exact RNA stanza from `.gitmodules`.
3. `git add` stages only that `.gitmodules` edit.
4. `git rm --cached` stages removal of only the exact `160000` gitlink; the checkout was already removed by `deinit`.
5. `rmdir` removes only an empty leftover mount directory and refuses a nonempty one.
6. The conditional removes the exact local `.git/config` section if this Git version left it behind.
7. `mv` removes the live `.git/modules/<path>` directory while retaining it as a reversible quarantine.

Do not substitute `git rm -f`, `rm -rf "$sub_path"`, or a hand edit of `.gitmodules` for these commands.

## Phase 3: verify and commit only the removal

```bash
test -z "$(git ls-files --stage -- "$sub_path")"
test -z "$(git config -f .gitmodules --get "submodule.${sub_name}.path")"
test ! -e "$sub_path"
test ! -e "$module_gitdir"
test -d "$retired_gitdir"
sciagent_after=$(git config -f .gitmodules --get-regexp '^submodule\.01_modules/SciAgent-toolkit\.')
test "$sciagent_after" = "$sciagent_before"
git diff --cached --check
git diff --cached --name-status
git diff --cached -- .gitmodules "$sub_path"
git status --short
```

Effects: the tests prove the gitlink, stanza, working tree, active administrative repository, and local config are gone; the quarantine remains; SciAgent's parsed config values equal the captured values; and the final commands expose the staged patch and complete dirty state for review. `git diff --cached --name-status` must name exactly `.gitmodules` and `$sub_path`. Compare `git status --short` with the Phase 1 baseline; all unrelated entries must be unchanged.

Also run:

```bash
git config --local --get-regexp "^submodule\.${sub_name}\."
rg -n --hidden --glob '!.git/**' --glob "!${retired_gitdir}/**" 'RNAseq-toolkit|01_(modules|scripts)/RNAseq-toolkit' .
```

Effect: the config query must have no output and exit 1. The search may find the deliberately staged deletion in `.gitmodules`; it must find no live outside reference or checkout.

Commit with an exact pathspec:

```bash
git commit --only -m 'Retire vendored RNAseq-toolkit submodule' -- .gitmodules "$sub_path"
removal_commit=$(git rev-parse HEAD)
git show --stat --oneline "$removal_commit"
git status --short
```

Effect: `--only` constructs the commit from the two named paths and cannot include unrelated dirty paths. The remaining commands record and review the new commit and worktree. Stop if the commit contains anything except the exact stanza and gitlink deletion.

Publish normally only after owner review:

```bash
repo_branch=$(git branch --show-current)
test -n "$repo_branch"
git push origin "$repo_branch"
```

Effect: this pushes the named current branch without rewriting history. Repository protections may require a feature branch and pull request instead; follow that repository's policy and do not force-push.

## Exact per-repository invocations

For each block, set the values, then run Phases 1–3. M1/M2/M3 and D1/D2/D3 remain subject to the owner hold.

M1:

```bash
repo_root='/data1/users/antonz/projects/Meta-Aging/14616-DM'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

M2:

```bash
repo_root='/data1/users/antonz/projects/Meta-Aging/14761-DM'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

M3:

```bash
repo_root='/data1/users/antonz/projects/Meta-Aging/14782-DM'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

M4:

```bash
repo_root='/data1/users/antonz/projects/Meta-Aging/neuroimmune-receptor-atlas'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

D1:

```bash
repo_root='/scratch/current/antonz/projects/DC-nexus/DC_Dictionary'
sub_path='01_scripts/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

D2:

```bash
repo_root='/scratch/current/antonz/projects/DC-nexus/DC_hum_verse'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

D3:

```bash
repo_root='/scratch/current/antonz/projects/DC-nexus/DC_mouse_cancer'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

After the child pointer updates described in the next section, run the same transaction for M5:

```bash
repo_root='/data1/users/antonz/projects/Meta-Aging'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

Then run it for D4:

```bash
repo_root='/scratch/current/antonz/projects/DC-nexus'
sub_path='01_modules/RNAseq-toolkit'
sub_name="$sub_path"
cd "$repo_root"
```

## Required parent gitlink updates

Removing the nested declaration creates a new commit in each child repository. A parent does not acquire that content change automatically: its gitlink must be advanced to the reviewed child commit. Do this only after every child removal commit is published and every child checkout is clean.

For `Meta-Aging`, after owners have resolved the three pre-existing pointer mismatches:

```bash
cd /data1/users/antonz/projects/Meta-Aging
git diff --cached --quiet
git -C 14616-DM status --porcelain=v1
git -C 14761-DM status --porcelain=v1
git -C 14782-DM status --porcelain=v1
git -C neuroimmune-receptor-atlas status --porcelain=v1
git add -- 14616-DM 14761-DM 14782-DM neuroimmune-receptor-atlas
git diff --cached --submodule=log -- 14616-DM 14761-DM 14782-DM neuroimmune-receptor-atlas
git diff --cached --name-status
git commit --only -m 'Advance subrepositories after RNAseq-toolkit retirement' -- 14616-DM 14761-DM 14782-DM neuroimmune-receptor-atlas
parent_branch=$(git branch --show-current)
test -n "$parent_branch"
git push origin "$parent_branch"
```

Effects: the first command requires an initially empty index; the four status commands must print nothing; `git add` stages only the four gitlink SHAs; the diffs let the owners approve the complete old-to-new commit ranges; `commit --only` records only those pointers; and the final commands publish the named parent branch without rewriting it. If a child did not change, omit it only after confirming its removal commit is already the recorded gitlink.

For `DC-nexus`, after owners have resolved all three pre-existing pointer mismatches:

```bash
cd /scratch/current/antonz/projects/DC-nexus
git diff --cached --quiet
git -C DC_Dictionary status --porcelain=v1
git -C DC_hum_verse status --porcelain=v1
git -C DC_mouse_cancer status --porcelain=v1
git add -- DC_Dictionary DC_hum_verse DC_mouse_cancer
git diff --cached --submodule=log -- DC_Dictionary DC_hum_verse DC_mouse_cancer
git diff --cached --name-status
git commit --only -m 'Advance subrepositories after RNAseq-toolkit retirement' -- DC_Dictionary DC_hum_verse DC_mouse_cancer
parent_branch=$(git branch --show-current)
test -n "$parent_branch"
git push origin "$parent_branch"
```

Effects are the same for the three DC child pointers. The unrelated `.devcontainer`, `.vscode`, `.claude`, `integration`, and other measured DC root changes remain unstaged and uncommitted. Compare the before/after `git status --short` output to prove that.

Review and publish each parent pointer commit normally. Only then run M5 or D4 respectively, so the root removal has an empty index and its own two-path commit.

## Phase 4: permanently remove the quarantined administrative repository

Keep each quarantine through review, publication, and any required pull-request merge. Before deleting it, re-run the remote-containment check for the captured `checkout_sha`, the local-ref and unreachable-object inventory, and confirm the removal commit is on the accepted remote branch.

**STOP POINT B — last lossless local rollback point:** do not continue unless the owner accepts that the exact local refs, reflogs, and objects inside the quarantine are no longer needed.

```bash
test -n "$retired_gitdir"
case "$retired_gitdir" in */modules/*RNAseq-toolkit.retired-rnaseq) true ;; *) false ;; esac
test -d "$retired_gitdir"
find "$retired_gitdir" -maxdepth 1 -mindepth 1 -print
rm -r -- "$retired_gitdir"
test ! -e "$retired_gitdir"
```

Effects: the first three commands validate the exact, narrow quarantine target; `find` displays its top level for final review; `rm -r` permanently deletes only that administrative repository; the final test proves it is gone. This is the first step that prevents an exact, offline restoration of local-only refs, reflogs, or objects.

## Rollback map

| Last completed step | Rollback |
|---|---|
| Phase 1 only | No worktree or index rollback is needed. `fetch` only refreshed remote-tracking refs. |
| `git submodule deinit` | While the stanza and gitlink remain, run `git submodule init -- "$sub_path"` and `git submodule update -- "$sub_path"`; reapply any captured custom local config keys. |
| `.gitmodules` stanza removed/staged | Because the gate required `.gitmodules` to be clean, run `git restore --staged --worktree -- .gitmodules` to restore it exactly. Confirm the SciAgent stanza still matches the capture. |
| Gitlink staged for deletion | Run `git restore --staged -- "$sub_path"`, restore `.gitmodules` as above, then initialise/update the checkout. |
| Local config section removed | Restore `.gitmodules` and the gitlink, run `git submodule init -- "$sub_path"`, then reapply any captured nonstandard `submodule.<name>.*` keys with `git config --local <key> <value>`. |
| Empty mount directory removed | Restoring and updating the submodule recreates it. |
| Administrative gitdir quarantined | Before initialising, run `mv -- "$retired_gitdir" "$module_gitdir"`; then restore the stanza/gitlink and run `git submodule update --init -- "$sub_path"`. |
| Removal staged but not committed | Move the quarantine back, run `git restore --staged -- "$sub_path" .gitmodules`, `git restore --worktree -- .gitmodules`, and `git submodule update --init -- "$sub_path"`; restore captured custom local config. |
| Removal committed, quarantine retained | Use `git revert <removal-commit>` rather than reset/amend, move the quarantine back, and run `git submodule update --init -- "$sub_path"`. Review before pushing the revert. |
| Parent gitlink update committed | Revert that exact parent commit; do not rewrite it. This returns the parent to the earlier child SHAs while retaining the child cleanup commits. |
| Quarantine permanently deleted | Revert the removal commit if desired, then initialise from the still-existing remote. This reconstructs tracked history, but exact local-only refs/reflogs/objects cannot be restored from the deleted directory. |

Publishing a removal commit is the first point after which returning a **public branch's commit sequence** to its former shape would require forbidden rewriting; the permitted rollback is an additive `git revert`. Permanently deleting the quarantine is the first point after which rollback is not an exact, lossless restoration of the local repository metadata.

## Historical record

Removing a checkout, gitlink, configuration stanza, and local `.git/modules` cache does not delete the independent `RNAseq-toolkit` repository or any commit already reachable there. Earlier superproject commits retain their historical gitlinks and can still be inspected; checking them out and initialising the old submodule continues to depend on the toolkit remote remaining available. The remote-containment gate prevents deletion when the checked-out commit is unique locally, while the separate local-ref/object review protects other local-only data.

This procedure adds ordinary removal and pointer-update commits only. It requires no rewriting of either superproject, any child repository, or the toolkit's published history. Rewriting published history is outside scope and forbidden without separately named consent.
