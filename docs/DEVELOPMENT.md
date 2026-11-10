# Development and review

Use the clean environment and commands in the [README](../README.md). The local reference
interpreter is the conda `20220806-03` Python 3.10.5, used only to create this repository's
isolated `.venv`. CI uses Python 3.10, 3.11 and 3.12.

## Milestone lifecycle

Invoke [the milestone skill](../.agents/skills/milestone/SKILL.md). It selects the next
eligible milestone, creates a `codex/<milestone>-<slug>` branch, implements and verifies
the work, opens a draft PR from the template, waits for final-revision Tests, and marks the
PR ready. The `ready_for_review` event requests one adversarial Claude pass. The agent waits
for the review, addresses findings, and requests another pass after material fixes by
repeating draft → green Tests → ready. Finished PRs stay unmerged until the researcher
instructs otherwise.

`Tests` runs on pushes and pull requests. Its 10- and 12-minute job timeouts cover
installation; the test processes themselves have separately enforced 120 s and 300 s
budgets (`scripts/run_tests.py`). Neither the workflow limit nor moving a test to the slow
tier relaxes the scientific test budget.

The Claude workflow triggers only on `ready_for_review` and `reopened`, so draft creation
and pushes do not spend a review. Same-repository non-draft PRs are eligible; fork PRs are
not reviewed with repository credentials.

## One-time Claude setup

The workflow file is installed, but a review cannot run until the repository has its own
credential; secrets are not copied between repositories.

1. Obtain a Claude Code OAuth token with `claude setup-token` on a machine signed in to the
   intended Claude account.
2. Add it as the repository secret `CLAUDE_CODE_OAUTH_TOKEN`, either under
   [Actions secrets](https://github.com/landreman/alpha_analysis_3/settings/secrets/actions)
   or with `gh secret set CLAUDE_CODE_OAUTH_TOKEN --repo landreman/alpha_analysis_3` at the
   interactive prompt. Never put the token in a file.
3. Ensure the [Claude GitHub App](https://github.com/apps/claude) has access to
   `landreman/alpha_analysis_3` (an existing installation may restrict repositories).
4. Open a real change as a draft PR, wait for Tests, then mark it ready, and confirm that
   `Claude Code Review` posts the top-level review comment.

The workflow fails with these instructions when the secret is absent; it never reports a
skipped reviewer as a successful review.

## Codex project

Use the repository folder as the project's primary folder so that `AGENTS.md`,
`.codex/config.toml` and `.agents/skills/milestone` are discovered. The project is
independent of the earlier repositories.

## Contents and provenance

The five boozmn files and the original note are byte-preserved and hashed
(`data/manifest.json`, `docs/references/provenance.json`). The Boozer reader was copied
from the original repository and is maintained here. The design-phase prototype lives in
`prototype/` with its evidence in `docs/validation/prototype/`. No large reference data,
mesh pipeline, certificate machinery or old work queue was copied; `docs/LESSONS.md`
records what the earlier attempts taught.
