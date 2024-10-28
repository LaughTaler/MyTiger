#!/bin/bash

# 设置 HTTP 账户名和密码的环境变量
#export HTTP_USERNAME="your_http_username"
#export HTTP_PASSWORD="your_http_password"

# 执行 expect 脚本
expect << EOF
spawn ./linux_build.sh

expect "Username:"
send "\$env(HTTP_USERNAME)\r"

expect "Password:"
send "\$env(HTTP_PASSWORD)\r"

interact
EOF
