# ����ԴĿ¼��Ŀ��Ŀ¼
$source_dir = "..\..\build\Debug"
$target_dir_lib = "..\..\dabao\meshrepair\meshrepair\windows\x86_64\lib"
$target_dir_bin = "..\..\dabao\meshrepair\meshrepair\windows\x86_64\bin"
$source_conanfile_dir = "..\..\dabao\meshrepair"

# ���Ŀ��Ŀ¼�����ڣ��򴴽�Ŀ¼
if (-not (Test-Path -Path $target_dir_lib)) {
    New-Item -ItemType Directory -Path $target_dir_lib
}

if (-not (Test-Path -Path $target_dir_bin)) {
    New-Item -ItemType Directory -Path $target_dir_bin
}

# ���� .lib �� .dll �ļ���Ŀ��Ŀ¼
Copy-Item -Path "$source_dir\*.lib" -Destination $target_dir_lib
Copy-Item -Path "$source_dir\*.dll" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.exe" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.pdb" -Destination $target_dir_bin
Copy-Item -Path "$source_dir\*.exp" -Destination $target_dir_bin

Write-Output "��ɸ����ļ�."


# ���� conanfile.py �� source_conanfile_dir
Copy-Item -Path ".\conanfile.py" -Destination $source_conanfile_dir

Set-Location $source_conanfile_dir

# ִ�� conan export-pkg ����
Invoke-Expression "conan export-pkg . -s os=`"Windows`" -s arch=`"x86_64`""
Invoke-Expression "conan upload meshrepair/1.0.0 -r=lab"
Pause
