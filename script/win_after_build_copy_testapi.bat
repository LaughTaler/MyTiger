@echo off
REM 复制文件
@REM copy build\release\Release\*.dll C:\Users\lab\tiger\testapi\
@REM copy build\release\Release\test_api.exe C:\Users\lab\tiger\testapi\

@REM REM 检查复制是否成功
@REM if %errorlevel% neq 0 (
@REM     echo 文件复制失败
@REM     exit /b %errorlevel%
@REM )

REM 切换到目标目录
cd /d C:\Users\lab\tiger\aft_tri_test\

REM 运行复制后的文件
AFT_Tri_test.exe [test1]
REM 检查执行是否成功
if %errorlevel% neq 0 (
    echo %errorlevel%
    echo failed to run exe
    exit /b %errorlevel%
)
AFT_Tri_test.exe [test2]
REM 检查执行是否成功
if %errorlevel% neq 0 (
    echo %errorlevel%
    echo failed to run exe
    exit /b %errorlevel%
)
AFT_Tri_test.exe [test3]
REM 检查执行是否成功
if %errorlevel% neq 0 (
    echo %errorlevel%
    echo failed to run exe
    exit /b %errorlevel%
)
AFT_Tri_test.exe [test4]
REM 检查执行是否成功
if %errorlevel% neq 0 (
    echo %errorlevel%
    echo failed to run exe
    exit /b %errorlevel%
)

@REM cd ..\dt_remesh_test\

@REM DT_Remesh_test.exe [test]
@REM REM 检查执行是否成功
@REM if %errorlevel% neq 0 (
@REM     echo %errorlevel%
@REM     echo failed to run exe
@REM     exit /b %errorlevel%
@REM )

echo operation success
