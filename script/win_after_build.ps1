# 定义源目录和目标目录
$SRC_DIR = "C:/Users/lab/tiger"
$DEST_DIR = "C:/Users/lab/tiger/libtiger/libtiger/windows/x86_64"
$CONANFILE_DEST_DIR = "C:/Users/lab/tiger/libtiger"

# 要复制的文件夹列表
$FOLDERS = @("lib", "bin", "include")

# 创建目标目录（如果不存在的话）
if (-Not (Test-Path -Path $DEST_DIR)) {
    New-Item -ItemType Directory -Force -Path $DEST_DIR
}

# 循环复制每个文件夹
foreach ($FOLDER in $FOLDERS) {
    $SRC_FOLDER = Join-Path -Path $SRC_DIR -ChildPath $FOLDER
    if (Test-Path -Path $SRC_FOLDER) {
        Copy-Item -Path $SRC_FOLDER -Destination $DEST_DIR -Recurse -Force
        Write-Output "Copied $FOLDER to $DEST_DIR"
    }
    else {
        Write-Output "Source folder $SRC_FOLDER does not exist"
    }
}

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

$command3 = "conan remove libtiger/$version -c -r=lab"

# 执行指令
Write-Output "执行指令: $command3"
Invoke-Expression $command3

$command4 = "conan upload libtiger/$version -r=lab"

# 执行指令
Write-Output "执行指令: $command4"
Invoke-Expression $command4
