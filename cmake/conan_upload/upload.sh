#!/bin/bash


# 定义源目录和目标目录
SRC_DIR=~/tiger
DEST_DIR=~/tiger/libtiger/libtiger/windows/x86_64
CONANFILE_DEST_DIR=~/tiger/libtiger

# 要复制的文件夹列表
FOLDERS=("lib" "tiger" "include")

# 创建目标目录（如果不存在的话）
mkdir -p "$DEST_DIR"

# 循环复制每个文件夹
for FOLDER in "${FOLDERS[@]}"; do
    # 检查源文件夹是否存在
    if [ -d "$SRC_DIR/$FOLDER" ]; then
        cp -r "$SRC_DIR/$FOLDER" "$DEST_DIR"
        echo "Copied $FOLDER to $DEST_DIR"
    else
        echo "Source folder $SRC_DIR/$FOLDER does not exist"
    fi
done

# 复制 conanfile.py 文件
CONANFILE_SRC="conanfile.py"
if [ -f "$CONANFILE_SRC" ]; then
    cp "$CONANFILE_SRC" "$CONANFILE_DEST_DIR"
    echo "Copied conanfile.py to $CONANFILE_DEST_DIR"
else
    echo "conanfile.py does not exist in $SRC_DIR"
fi



cd "$CONANFILE_DEST_DIR" || { echo "Failed to change directory to $CONANFILE_DEST_DIR"; exit 1; }
echo "Changed directory to $CONANFILE_DEST_DIR"


echo "All done!"













# 获取当前月份和日期
current_month=$(date +%-m)  # 2位数字的月份，如“01”表示一月
current_day=$(date +%-d)    # 2位数字的日期，如“01”表示一号

# 构建用于conan remove的版本号
version="${current_month}.${current_day}"

# 构建完整的指令
command="conan remove libtiger/${version} -c"

# 执行指令
echo "执行指令: ${command}"
${command}


command2='conan export-pkg . -s os='Linux' -s arch='x86_64''

# 执行指令
echo "执行指令: ${command2}"
${command2}

command3="conan remove libtiger/${version} -c -r=lab"

# 执行指令
echo "执行指令: ${command3}"
${command3}

command4="conan upload libtiger/${version} -r=lab"

# 执行指令
echo "执行指令: ${command4}"
${command4}