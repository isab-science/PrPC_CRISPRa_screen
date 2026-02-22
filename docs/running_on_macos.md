# Instruction for running on a Macintosh

This guide assumes you are new to macOS command line usage.

Use the project via the web interface (`webapp`).  
Do not use `orchestrate_screen_workflow.ps1` on macOS (that GUI is Windows-only).

## 1) Open the Terminal app

1. Press `Command + Space`.
2. Type `Terminal`.
3. Press `Enter`.

You will paste all commands below into this Terminal window.

## 2) Install Apple Command Line Tools

In Terminal, run:

```bash
xcode-select --install
```

A popup may appear. Click Install and wait for completion.

## 3) Install Homebrew (only if missing)

Check whether Homebrew is already installed:

```bash
brew --version
```

- If this prints a version number, continue to step 4.
- If it says command not found, install Homebrew with:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

After install, close Terminal, open Terminal again, and run:

```bash
brew --version
```

## 4) Install Python

In Terminal:

```bash
brew install python@3.11
python3 --version
```

You should see something like `Python 3.11.x`.

## 5) Download the project

In Terminal, move to your home folder, clone the repository, and enter it:

```bash
cd ~
git clone https://github.com/isab-science/PrPC_CRISPRa_screen.git
cd PrPC_CRISPRa_screen
```

## 6) Create a project virtual environment

In Terminal:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

When activated, your prompt usually starts with `(.venv)`.

## 7) Install required Python packages

In the same Terminal window:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 8) Start the web app

In Terminal:

```bash
python -m uvicorn webapp.app:app --reload
```

Leave this Terminal window open while using the app.

## 9) Open the app in your browser

Open Safari (or Chrome) and go to:

```text
http://127.0.0.1:8000
```

## 10) Run the pipeline from the web page

1. Set `Data root` to the folder that contains your screen files.
2. Choose `Raw dir/file`.
3. Choose `Layout CSV`.
4. Choose `Genomics XLSX`.
5. Choose `Output dir`.
6. Review optional fields (`Sheet`, `Skip FRET`, `Skip GLO`, `Heatmap plate`).
7. Click `Run Pipeline`.

Generated outputs are written into your selected output directory.

## 11) Next time you run it

Only run:

```bash
cd ~/PrPC_CRISPRa_screen
source .venv/bin/activate
python -m uvicorn webapp.app:app --reload
```

## Common issues

- `ModuleNotFoundError`:
  Run `python -m pip install -r requirements.txt` again after activating `.venv`.
- Port already in use:
  Start on another port:
  `python -m uvicorn webapp.app:app --reload --port 8001`
  then open `http://127.0.0.1:8001`.
- Wrong window:
  All commands must be run in the macOS `Terminal` app, not in Finder, not in Word, not in browser.
