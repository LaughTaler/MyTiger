# 定义源目录和目标目录
$source_dir = "..\..\build\Debug"
$target_dir_lib = "..\..\dabao\meshrepair\meshrepair\windows\x86_64\lib"
$target_dir_bin = "..\..\dabao\meshrepair\meshrepair\windows\x86_64\bin"
$source_conanfile_dir = "..\..\dabao\meshrepair"

# 如果目标目录不存在，则创建目录
if (-not (Test-Path -Path $target_dir_lib)) {
    New-Item -ItemType Directory -Path $target_dir_lib
}

if (-not (Test-Path -Path $target_dir_bin)) {
    New-Item -ItemType Directory -Path $target_dir_bin
}

# 复制 .lib 和 .dll 文件到目标目录
Copy-Item -Path "$source_dir\*.lib" -Destination $target_dir_lib
Copy-Item -Path "$source_dir\*.dll" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.exe" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.pdb" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.exp" -Destination $target_dir_bin

Write-Output "完成复制文件."


# 复制 conanfile.py 到 source_conanfile_dir
Copy-Item -Path ".\conanfile.py" -Destination $source_conanfile_dir

Set-Location $source_conanfile_dir

# 执行 conan export-pkg 命令
Invoke-Expression "conan export-pkg . -s os=`"Windows`" -s arch=`"x86_64`""
Invoke-Expression "conan upload meshrepair/1.0.0 -r=lab"
Pause
