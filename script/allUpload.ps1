# 定义源目录和目标目录
$SRC_DIR = "C:/Users/lab/tiger"
$DEST_DIR = "C:/Users/lab/tiger/libtiger/libtiger/windows/x86_64"
$CONANFILE_DEST_DIR = "C:/Users/lab/tiger/libtiger"
$DEST_DIR2 = "C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64"

$DEST_DIR1d = "C:/Users/lab/tiger/libtiger/libtiger/windows/x86_64/all/debug"
$DEST_DIR1r = "C:/Users/lab/tiger/libtiger/libtiger/windows/x86_64/all/release"
$DEST_DIR2d = "C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/debug"
$DEST_DIR2r = "C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/release"
$SRC_DIRd = "C:/Users/lab/tiger/debug"
$SRC_DIRr = "C:/Users/lab/tiger/release"
# 要复制的文件夹列表
$FOLDERS = @("lib", "bin", "include")

if (Test-Path -Path $DEST_DIR) {
    # 如果目录存在，删除该目录
    Remove-Item -Path $DEST_DIR -Recurse -Force
}

# 创建目标目录
New-Item -ItemType Directory -Force -Path $DEST_DIR
New-Item -ItemType Directory -Force -Path $DEST_DIR1d
New-Item -ItemType Directory -Force -Path $DEST_DIR1r

if (Test-Path -Path $DEST_DIR2) {
    # 如果目录存在，删除该目录
    Remove-Item -Path $DEST_DIR2 -Recurse -Force
}

# 创建目标目录
New-Item -ItemType Directory -Force -Path $DEST_DIR2
New-Item -ItemType Directory -Force -Path $DEST_DIR2d
New-Item -ItemType Directory -Force -Path $DEST_DIR2r
# 循环复制每个文件夹
foreach ($FOLDER in $FOLDERS) {
    $SRC_FOLDERd = Join-Path -Path $SRC_DIRd -ChildPath $FOLDER
    if (Test-Path -Path $SRC_FOLDERd) {
        Copy-Item -Path $SRC_FOLDERd -Destination $DEST_DIR1d -Recurse -Force
        Write-Output "Copied $FOLDER to $DEST_DIR1d"
    }
    else {
        Write-Output "Source folder $SRC_FOLDERd does not exist"
    }
}
# 循环复制每个文件夹
foreach ($FOLDER in $FOLDERS) {
    $SRC_FOLDERr = Join-Path -Path $SRC_DIRr -ChildPath $FOLDER
    if (Test-Path -Path $SRC_FOLDERr) {
        Copy-Item -Path $SRC_FOLDERr -Destination $DEST_DIR1r -Recurse -Force
        Write-Output "Copied $FOLDER to $DEST_DIR1r"
    }
    else {
        Write-Output "Source folder $SRC_FOLDERr does not exist"
    }
}

scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/debug/bin /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/debug/
scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/debug/include /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/debug/
scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/debug/lib /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/debug/
scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/release/bin /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/release/
scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/release/include /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/release/
scp -i C:\Users\lab\.ssh\id_rsa -r user@192.168.0.3:~/../gitlab-runner/tiger/release/lib /C:/Users/lab/tiger/libtiger/libtiger/linux/x86_64/all/release/


# # 复制 conanfile.py 文件
# $CONANFILE_SRC = "conanfile.py"
# if (Test-Path -Path $CONANFILE_SRC) {
#     Copy-Item -Path $CONANFILE_SRC -Destination $CONANFILE_DEST_DIR -Force
#     Write-Output "Copied conanfile.py to $CONANFILE_DEST_DIR"
# } else {
#     Write-Output "conanfile.py does not exist in $CONANFILE_SRC"
# }

Set-Location -Path $CONANFILE_DEST_DIR -ErrorAction Stop
Write-Output "Changed directory to $CONANFILE_DEST_DIR"

Write-Output "All done!"











# 获取当前月份和日期
$current_month = (Get-Date).Month
$current_day = (Get-Date).Day

# 构建用于conan remove的版本号
$version = "$current_month.$current_day"

# 构建完整的指令
$command = "conan remove libtiger/$version -c"

# 执行指令
Write-Output "执行指令: $command"
Invoke-Expression $command

$command2 = "conan export-pkg . -s os='Windows' -s arch='x86_64'"

# 执行指令
Write-Output "执行指令: $command2"
Invoke-Expression $command2

$command22 = "conan export-pkg . -s os='Linux' -s arch='x86_64'"

# 执行指令
Write-Output "执行指令: $command22"
Invoke-Expression $command22

$command3 = "conan remove libtiger/$version -c -r=lab"

# 执行指令
Write-Output "执行指令: $command3"
Invoke-Expression $command3

$command4 = "conan upload libtiger/$version -r=lab"

# 执行指令
Write-Output "执行指令: $command4"
Invoke-Expression $command4
