@echo off
setlocal

REM ����ԴĿ¼��Ŀ��Ŀ¼
set "source_dir=..\..\build\Debug"
set "target_dir_lib=..\..\dabao\meshrepair\meshrepair\windows\x86_64\lib"
set "target_dir_bin=..\..\dabao\meshrepair\meshrepair\windows\x86_64\bin"
set "source_conanfile_dir=..\..\dabao\meshrepair"
REM ���Ŀ��Ŀ¼�����ڣ��򴴽�Ŀ¼
if not exist "%target_dir_lib%" (
    mkdir "%target_dir_lib%"
)
if not exist "%target_dir_bin%" (
    mkdir "%target_dir_bin%"
)

REM ���� .lib �� .dll �ļ���Ŀ��Ŀ¼
copy "%source_dir%\*.lib" "%target_dir_lib%"
copy "%source_dir%\*.dll" "%target_dir_bin%"
copy "%source_dir%\*.exe" "%target_dir_bin%"
copy "%source_dir%\*.pdb" "%target_dir_bin%"
copy "%source_dir%\*.exp" "%target_dir_bin%"

echo ��ɸ����ļ�.
endlocal
pause
