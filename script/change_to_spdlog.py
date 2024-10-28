import os
import re
import sys

def change_printf(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    replaced_content = content
    # 匹配 printf 语句的正则表达式
    printf_pattern = r'printf\("([^"]+)"(, .*)?\);'

    # 替换成 spdlog::info 的模板
    spdlog_template = 'spdlog::info("{}", {});'

    # 匹配 printf 语句并进行替换
    result = re.sub(printf_pattern, lambda m: spdlog_template.format(m.group(1).replace('%d', '{}').replace('%f', '{}'), m.group(2) if m.group(2) else ''), replaced_content)



def add_spdlog_include(file_path):
    # 读取文件内容
    with open(file_path, 'r') as file:
        content = file.readlines()

    # 找到第一个非空行的位置
    include_index = 0
    for i, line in enumerate(content):
        if line.strip().startswith("#include \""):
            include_index = i
            break
    if include_index == 0:
        for i, line in enumerate(content):
            if line.strip().startswith("#include <"):
                include_index = i
                break
    # 在第一个非空行之后插入spdlog的include
    content.insert(include_index, "#include <spdlog/spdlog.h> \n ")

    #print(include_index)
    # 将内容写回文件
    #with open(file_path, 'w') as file:
    #    file.writelines(content)

def replace_cout_and_printf_with_spdlog_in_file(file_path):
    # 读取文件内容
    with open(file_path, 'r') as file:
        content = file.read()
    replaced_content = content
    # 使用正则表达式进行替换
    #replaced_content = re.sub(r'std::cout\s*<<\s*"([^"<]+)"\s*<<\s*std::endl\s*;', r'spdlog::info("\1")', replaced_content)
    #replaced_content = re.sub(r'std::cout\s*<<\s*"([^"]+)"\s*;', r'spdlog::info("\1")', replaced_content)
    #replaced_content = re.sub(r'printf\s*\(\s*"([^"]+)"\s*\);', r'spdlog::info("\1")', replaced_content)

    replaced_content = re.sub(r'std::cout\s*<<\s*([^"]+)\s*;', r'spdlog::info("\1")', replaced_content)
    replaced_content = re.sub(r'std::cout\s*\<\<\s*([^\"]+)\s*\<\<\s*std::endl\s*;', r'spdlog::info("\1")', replaced_content)
    replaced_content = re.sub(r'std::cout\s*<<\s*(.*?)\s*<<\s*std::endl;', r'spdlog::info(\1);', replaced_content)
    replaced_content = re.sub(r'std::cout\s*<<\s*(.*?)<<\s*"\\n"\s*;', r'spdlog::info(\1);', replaced_content)
    replaced_content = re.sub(r'std::cout\s*<<\s*(.*?)\\n\s*;', r'spdlog::info(\1);', replaced_content)
    replaced_content = re.sub(r'std::cout\s*<<\s*(.*?)\s*;', r'spdlog::info(\1);', replaced_content)
    replaced_content = re.sub(r'printf\s*\(\s*([^"]+)\s*\);', r'spdlog::info("\1")', replaced_content)

     # 匹配 printf 语句的正则表达式
    printf_pattern = r'printf\("([^,"]+)"(, .*)?\);'

    # 替换成 spdlog::info 的模板
    spdlog_template = 'spdlog::info("{}", {});'

    # 匹配 printf 语句并进行替换
    replaced_content = re.sub(printf_pattern, lambda m: spdlog_template.format(m.group(1).replace('%d', '{}').replace('%f', '{}'), m.group(2) if m.group(2) else ''), replaced_content)

    # 匹配 spdlog 语句中的额外逗号的正则表达式
    extra_comma_pattern = r'(spdlog::\w+)\("([^"]+)"\s*,\s*\);'

    # 去除额外逗号的模板
    fixed_spdlog_template = r'\1("\2");'

    # 匹配 spdlog 语句并去除额外逗号
    replaced_content = re.sub(extra_comma_pattern, fixed_spdlog_template, replaced_content)

    # 匹配含有多余逗号的 spdlog::info 语句的正则表达式
    extra_comma_pattern = r'(spdlog::\w+)\("([^"]+)",\s*,(.+)\);'

    # 去除多余逗号的模板
    fixed_spdlog_template = r'\1("\2",\3);'

    # 匹配 spdlog::info 语句并去除多余逗号
    replaced_content = re.sub(extra_comma_pattern, fixed_spdlog_template, replaced_content)


    #print(replaced_content)
    # 将替换后的内容写回文件
    with open(file_path, 'w') as file:
        file.write(replaced_content)

def process_file(file_path):
    # 只处理C++源代码文件
    if file_path.endswith(".cpp")  or file_path.endswith(".hpp"):
        print("prossing",file_path)
        add_spdlog_include(file_path)
        replace_cout_and_printf_with_spdlog_in_file(file_path)

def process_directory(directory):
    # 遍历目录下所有文件和子目录
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            process_file(file_path)

if __name__ == "__main__":

    
     # 获取命令行参数
    arguments = sys.argv

    directory_path = arguments[1]

    process_directory(directory_path)
