# Omniculator installer for Windows.
#   irm https://raw.githubusercontent.com/turtle170/Omniculator/main/install.ps1 | iex
# Downloads the latest release into %LOCALAPPDATA%\Programs\Omniculator and
# adds that folder to your user PATH.

& {
    $ErrorActionPreference = 'Stop'
    $repo = 'turtle170/Omniculator'
    $dir = Join-Path $env:LOCALAPPDATA 'Programs\Omniculator'
    $exe = Join-Path $dir 'omniculator.exe'

    Write-Host 'Installing Omniculator...'
    $release = Invoke-RestMethod "https://api.github.com/repos/$repo/releases/latest" `
        -Headers @{ 'User-Agent' = 'omniculator-installer' }
    $asset = $release.assets | Where-Object { $_.name -eq 'omniculator.exe' } | Select-Object -First 1
    if (-not $asset) { throw "The latest release ($($release.tag_name)) has no omniculator.exe." }

    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    Invoke-WebRequest $asset.browser_download_url -OutFile $exe -UseBasicParsing

    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    $parts = @($userPath -split ';' | Where-Object { $_ })
    if ($parts -notcontains $dir) {
        [Environment]::SetEnvironmentVariable('Path', (($parts + $dir) -join ';'), 'User')
        Write-Host "Added $dir to your PATH."
    }
    if (($env:Path -split ';') -notcontains $dir) { $env:Path += ";$dir" }

    Write-Host "Installed Omniculator $($release.tag_name) to $exe"
    Write-Host 'Try: omniculator "x^2 - 5x + 6 = 0"'
}
