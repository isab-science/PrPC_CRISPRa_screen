param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)

$ErrorActionPreference = 'Stop'

$hookFile = Join-Path $RepoRoot '.githooks\prepare-commit-msg'
if (-not (Test-Path $hookFile)) {
    throw "Hook file not found: $hookFile"
}

Push-Location $RepoRoot
try {
    git config core.hooksPath .githooks
    $configured = git config --get core.hooksPath
    if ($configured -ne '.githooks') {
        throw "Failed to configure core.hooksPath. Current value: $configured"
    }
    Write-Host 'Git hooks configured successfully.'
    Write-Host 'core.hooksPath=.githooks'
} finally {
    Pop-Location
}