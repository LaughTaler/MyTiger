@echo off
setlocal

REM 定义源目录和目标目录
set "source_dir=..\..\build\Debug"
set "target_dir_lib=..\..\dabao\meshrepair\meshrepair\windows\x86_64\lib"
set "target_dir_bin=..\..\dabao\meshrepair\meshrepair\windows\x86_64\bin"
set "source_conanfile_dir=..\..\dabao\meshrepair"
REM 如果目标目录不存在，则创建目录
if not exist "%target_dir_lib%" (
    mkdir "%target_dir_lib%"
)
if not exist "%target_dir_bin%" (
    mkdir "%target_dir_bin%"
)

REM 复制 .lib 和 .dll 文件到目标目录
copy "%source_dir%\*.lib" "%target_dir_lib%"
copy "%source_dir%\*.dll" "%target_dir_bin%"
copy "%source_dir%\*.exe" "%target_dir_bin%"
copy "%source_dir%\*.pdb" "%target_dir_bin%"
copy "%source_dir%\*.exp" "%target_dir_bin%"

echo 完成复制文件.
endlocal
pause
