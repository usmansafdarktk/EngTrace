# Runs one full-run stage under a Windows keep-awake request and appends its output to a log.
#
#   powershell -NoProfile -ExecutionPolicy Bypass -File full_run_28092026\keepawake_run.ps1 `
#       -Log full_run_28092026\scores\e5.log -m full_run_28092026.judge --yes --max-usd 55 --workers 16
#
# The request (ES_CONTINUOUS | ES_SYSTEM_REQUIRED) stops idle sleep for as long as this process lives.
# It changes no power setting and does not stop sleep from the lid, the power button or a flat battery,
# so the laptop still has to be on mains power with the lid open. The log gets a header line (time,
# process id, whether the request was granted, the command) and a footer line with the exit code.
param(
    [Parameter(Mandatory = $true)][string]$Log,
    [string]$Python = (Get-Command python -ErrorAction SilentlyContinue).Source,
    [Parameter(ValueFromRemainingArguments = $true)][string[]]$PyArgs
)
if (-not $Python) { throw 'python not found: pass -Python <path to python.exe>' }
Add-Type -Namespace KeepAwake -Name Native -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
$granted = [KeepAwake.Native]::SetThreadExecutionState([uint32]2147483649)
$env:PYTHONIOENCODING = 'utf-8'
$env:PYTHONUNBUFFERED = '1'
Set-Location (Split-Path $PSScriptRoot -Parent)
if (-not [IO.Path]::IsPathRooted($Log)) { $Log = [IO.Path]::GetFullPath((Join-Path (Get-Location) $Log)) }
$line = "`"$Python`" $($PyArgs -join ' ')"
$stamp = (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ')
"=== $stamp pid $PID keep-awake $(if ($granted) { 'held' } else { 'NOT GRANTED' }): $line" | Out-File -Append -Encoding ascii $Log
cmd.exe /c "`"$line >> `"$Log`" 2>&1`""
$code = $LASTEXITCODE
$stamp = (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ')
"=== $stamp exit $code" | Out-File -Append -Encoding ascii $Log
