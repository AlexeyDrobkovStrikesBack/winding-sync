# Publishing winding-sync — step by step

Written for your setup: GitHub user `abundantjoe`, Windows, PowerShell.

---

## The thing you asked earlier

**You don't paste files into a release.** Two separate things:

| | what it is | how files get there |
|---|---|---|
| **The repository** | your living code, browsable file by file | you *push* it with `git` |
| **A release** | a frozen, named snapshot of the repo | GitHub builds it **automatically** from a tag |

Push first, then tag a release pointing at it. GitHub generates the download
archives for you. You attach files only for things the repo doesn't contain — a
compiled binary, a trained model. You have neither, so attach nothing.

---

## Step 0 — Put the folder somewhere sensible

Extract `winding-sync.zip` to your Desktop. You should end up with:

```
C:\Users\Joseph Balmaceda\Desktop\winding-sync\
```

containing `winding_sync\`, `run_winding_test.py`, `README.md`, `LICENSE`,
`requirements.txt`, `.gitignore`, and `docs\img\winding_field.png`.

**This is not the Vesuvius folder.** Keep that one where it is — it holds the
triage tooling and your virtual environment, and it stays local.

## Step 1 — Install Git

Skip if `git --version` already prints something.

Download from https://git-scm.com/download/win, install with all defaults, then
**close and reopen PowerShell** (the PATH doesn't refresh in an open window):

```
git --version
```

## Step 2 — Identify yourself to Git

```
git config --global user.name "Joseph Balmaceda"
git config --global user.email "joseph.balmaceda@gmail.com"
```

**One privacy decision before you run that.** Whatever email you commit with
becomes permanently visible in the public commit log, and address scrapers do
read GitHub. If you'd rather not publish your personal address, GitHub gives
you a free alias: go to https://github.com/settings/emails, tick **Keep my
email addresses private**, and it shows an address like
`12345678+abundantjoe@users.noreply.github.com`. Then use:

```
git config --global user.email "YOUR-ID+abundantjoe@users.noreply.github.com"
```

Either choice works and commits still link to your account. Decide now —
changing it later doesn't rewrite commits you've already pushed.

## Step 3 — Create the empty repository

1. Go to https://github.com/new
2. **Repository name:** `winding-sync`
3. **Description:**
   ```
   Automatic relative winding constraints from CT for Herculaneum scroll unrolling, reconciled by L1 integer synchronization
   ```
4. Select **Public**
5. **Leave all three checkboxes unticked** — "Add a README file", "Add
   .gitignore", "Choose a license". You already have all three. Ticking any of
   them creates a conflicting commit you'd have to untangle on your first push.
6. Click **Create repository**

Your repo will live at **https://github.com/abundantjoe/winding-sync**

## Step 4 — Push the code

```
cd "C:\Users\Joseph Balmaceda\Desktop\winding-sync"
git init
git add .
git status
```

**Read the `git status` output before going further.** You should see about 15
files: the `winding_sync/` modules, `run_winding_test.py`, `README.md`,
`LICENSE`, `requirements.txt`, `docs/img/winding_field.png`.

You should **not** see `.venv` or `__pycache__`. If you do, stop and tell me —
`.venv` is hundreds of megabytes of machine-specific packages and must never go
into a repo.

Then:

```
git commit -m "v0.1.0: automatic winding constraints from CT via L1 synchronization"
git branch -M main
git remote add origin https://github.com/abundantjoe/winding-sync.git
git push -u origin main
```

A browser window opens to sign in to GitHub the first time. That's expected.

Refresh https://github.com/abundantjoe/winding-sync — your code is there and
the README renders with the banner image at the top.

## Step 5 — Create the release

1. Go to https://github.com/abundantjoe/winding-sync/releases/new
2. Click **Choose a tag**, type `v0.1.0`, then click
   **+ Create new tag: v0.1.0 on publish**
3. **Release title:**
   ```
   v0.1.0 — Automatic winding constraints from CT, reconciled by L1 synchronization
   ```
4. **Description:** open `RELEASE_NOTES.md` and paste everything under
   *"Release description (paste this into the GitHub release body)"*
5. **Attach nothing.** Leave the file-upload area empty.
6. Leave **Set as a pre-release** unticked — `0.1.0` already signals early
   status through the version number itself.
7. Click **Publish release**

GitHub now shows a release with auto-generated source archives.

## Step 6 — Submit to the Progress Prize

1. https://forms.gle/xoF5C3QsYutKP97x7
2. Paste the sections from `SUBMISSION_form_answers.md` into the matching
   fields
3. Repo URL: `https://github.com/abundantjoe/winding-sync`
4. **Deadline: 11:59pm Pacific, July 31st, 2026**

## Step 7 — Announce it

The Challenge explicitly weighs whether tools *get used*, and a repo nobody
knows about earns nothing. Post in the Discord (https://discord.gg/V4fJhvtaQn):
link the repo, two sentences on what it does, and say plainly that absolute
counts are uncalibrated and you'd welcome help. That reads as someone to
collaborate with rather than someone selling — and given the 0.66 ceiling
neither of us cracked, somebody there may see it quickly.

---

## Troubleshooting

**`git push` rejected — "remote contains work that you do not have locally"**
You ticked a checkbox in Step 3.
```
git pull origin main --allow-unrelated-histories
git push -u origin main
```

**`.venv` or `__pycache__` appearing in `git status`**
The `.gitignore` isn't in the folder, or Git tracked them before it existed.
```
git rm -r --cached .venv
git commit -m "Stop tracking local environment files"
```

**"fatal: not a git repository"**
You're in the wrong folder. Re-run the `cd` from Step 4.

**"support for password authentication was removed"**
An old Git version trying to use a password. Update Git, or install GitHub CLI
(https://cli.github.com) and run `gh auth login`.

**Banner image doesn't render on GitHub**
Confirm `docs/img/winding_field.png` actually pushed. GitHub paths are
case-sensitive: `docs/img/` is not `Docs/Img/`.

**Start over completely**
Delete the hidden `.git` folder inside `winding-sync` and redo Step 4. Nothing
outside that folder is affected.

---

## Quick reference

| | |
|---|---|
| Repo | https://github.com/abundantjoe/winding-sync |
| New release | https://github.com/abundantjoe/winding-sync/releases/new |
| Tag | `v0.1.0` |
| Submission form | https://forms.gle/xoF5C3QsYutKP97x7 |
| Deadline | 11:59pm Pacific, July 31st, 2026 |
| Discord | https://discord.gg/V4fJhvtaQn |
