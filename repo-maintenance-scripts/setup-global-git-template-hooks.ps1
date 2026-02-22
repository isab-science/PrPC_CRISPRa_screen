param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path,
    [string]$TemplateDir = (Join-Path $HOME '.git-template')
)

$ErrorActionPreference = 'Stop'

$sourceHook = Join-Path $RepoRoot '.githooks\prepare-commit-msg'
if (-not (Test-Path $sourceHook)) {
    throw "Source hook not found: $sourceHook"
}

$targetHooksDir = Join-Path $TemplateDir 'hooks'
$targetHook = Join-Path $targetHooksDir 'prepare-commit-msg'

if (-not (Test-Path $targetHooksDir)) {
    New-Item -ItemType Directory -Path $targetHooksDir -Force | Out-Null
}

Copy-Item -Path $sourceHook -Destination $targetHook -Force

# Ensure template hook is executable on POSIX environments too.
try {
    & git update-index --chmod=+x -- "$sourceHook" 2>$null | Out-Null
} catch {
}

git config --global init.templateDir "$TemplateDir"

$configured = git config --global --get init.templateDir
if ($configured -ne $TemplateDir) {
    throw "Failed to configure init.templateDir. Current value: $configured"
}

Write-Host 'Global git template hooks configured.'
Write-Host "init.templateDir=$TemplateDir"
Write-Host "Template hook: $targetHook"
Write-Host 'Applies automatically to future git init/clone operations.'