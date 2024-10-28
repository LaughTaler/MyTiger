@echo off

:: rmdir /s /q build -rf
mkdir build
cd build
cd release

whoami

:: 设置Visual Studio 2022环境变量
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

:: 运行CMake配置，指定生成器为Visual Studio 2022
cmake -G "Visual Studio 17 2022" -A x64 -DCMAKE_INSTALL_PREFIX=%USERPROFILE%/tiger/release ../.. -Wno-dev > compile_output.log

if %ERRORLEVEL% neq 0 (
    echo CMake配置失败
    exit /b %ERRORLEVEL%
)

:: 构建项目，使用多线程
cmake --build . --config Release  -- /m   >> compile_output.log

if %ERRORLEVEL% neq 0 (
    echo CMak config release失败
    exit /b %ERRORLEVEL%
)

:: 安装
cmake --install . --config Release  >> compile_output.log
if %ERRORLEVEL% neq 0 (
    echo CMake install release失败
    exit /b %ERRORLEVEL%
)

cd ../debug

:: 设置Visual Studio 2022环境变量
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

:: 运行CMake配置，指定生成器为Visual Studio 2022
cmake -G "Visual Studio 17 2022" -A x64 -DCMAKE_INSTALL_PREFIX=%USERPROFILE%/tiger/debug ../.. -Wno-dev > compile_output.log

if %ERRORLEVEL% neq 0 (
    echo CMake配置失败
    exit /b %ERRORLEVEL%
)

:: 构建项目，使用多线程
cmake --build . --config debug  -- /m  >> compile_output.log 
if %ERRORLEVEL% neq 0 (
    echo CMak build debug失败
    exit /b %ERRORLEVEL%
)

:: 安装
cmake --install . --config debug >> compile_output.log 
if %ERRORLEVEL% neq 0 (
    echo CMake install build失败
    exit /b %ERRORLEVEL%
)
cd ../..
if not exist "C:\Users\lab\tiger\aft_tri_test\" (
    mkdir "C:\Users\lab\tiger\aft_tri_test\"
)
copy build\release\Release\*.dll C:\Users\lab\tiger\aft_tri_test\
copy build\release\Release\AFT_Tri_test.exe C:\Users\lab\tiger\aft_tri_test\

@REM copy build\release\Release\*.dll C:\Users\lab\tiger\aft_tri_test\
@REM copy build\release\Release\AFT_Tri_test.exe C:\Users\lab\tiger\aft_tri_test\

@REM REM 检查复制是否成功
if %errorlevel% neq 0 (
    echo copy files failed
    exit /b %errorlevel%
)
@echo on
