"C:\Program Files\Git\bin\bash.exe" openblas-download.sh
powershell.exe -nologo -noprofile -command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('openblas-win32-int32.zip','.'); }"
powershell.exe -nologo -noprofile -command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('openblas-win32-int64.zip','.'); }"
